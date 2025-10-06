/* Stellarium Web Engine - Copyright (c) 2022 - Stellarium Labs SRL
 *
 * This program is licensed under the terms of the GNU AGPL v3, or
 * alternatively under a commercial licence.
 *
 * The terms of the AGPL v3 license can be found in the main directory of this
 * repository.
 */

#include "observer.h"
#include "swe.h"
#include "algos/utctt.h"
#include "navigation.h"

// Simple xor hash function.
static uint32_t hash_xor(uint32_t v, const char *data, int len)
{
    int i;
    // Only works on 4 bytes aligned data, either of size exactly equal to
    // 1 byte, or with a size that is a multiple of 4 bytes.
    assert((uintptr_t)data % 4 == 0);
    assert(len == 1 || len % 4 == 0);
    if (len == 1) return v ^ *data;
    for (i = 0; i < len; i += 4) {
        v ^= *(uint32_t*)(__builtin_assume_aligned(data, 4) + i);
    }
    return v;
}

static void update_matrices(observer_t *obs)
{
    eraASTROM *astrom = &obs->astrom;
    // We work with 3x3 matrices, so that we can use the erfa functions.
    double rdir[3][3];
    double ro2v[3][3];  // Rotate from observed to view.
    double ri2h[3][3];  // Equatorial J2000 (ICRF) to horizontal.
    double rh2i[3][3];  // Horizontal to Equatorial J2000 (ICRF).
    double ri2v[3][3];  // Equatorial J2000 (ICRF) to view.
    double ri2e[3][3];  // Equatorial J2000 (ICRF) to ecliptic.
    double re2i[3][3];  // Eclipic to Equatorial J2000 (ICRF).
    double rc2v[3][3];
    double view_rot[3][3];
    // r2gl changes the coordinate from z up to y up orthonomal.
    const double r2gl[3][3] = {{ 0, 0, -1},
                               {-1, 0,  0},
                               { 0, 1,  0}};

    const double flip_y[3][3] = {{1,  0, 0},
                                 {0, -1, 0},
                                 {0,  0, 1}};

    mat3_set_identity(rdir);
    mat3_rx(obs->roll, rdir, rdir);
    mat3_ry(obs->pitch, rdir, rdir);
    mat3_rz(-obs->yaw, rdir, rdir);

    if (mat3_det(obs->ro2m) > 0)
        mat3_product(ro2v, 4, r2gl, flip_y, rdir, obs->ro2m);
    else
        mat3_product(ro2v, 3, r2gl, rdir, obs->ro2m);

    // Extra rotation for screen center offset.
    assert(!isnan(obs->view_offset_alt));
    mat3_set_identity(view_rot);
    mat3_rx(obs->view_offset_alt, view_rot, view_rot);
    mat3_mul(view_rot, ro2v, ro2v);

    // Compute rotation matrix from CIRS to horizontal.
    mat3_set_identity(ri2h);
    // Earth rotation.
    mat3_rz(astrom->eral, ri2h, ri2h);
    // Polar motion.
    double rpl[3][3] = {{1, 0, astrom->xpl},
                        {0, 1, astrom->ypl},
                        {astrom->xpl, astrom->ypl, 1}};
    mat3_mul(ri2h, rpl, ri2h);
    // Cartesian -HA,Dec to Cartesian Az,El (S=0,E=90).
    mat3_ry(-obs->phi + M_PI / 2, ri2h, ri2h);
    double rsx[3][3] = {{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    mat3_mul(ri2h, rsx, ri2h);
    mat3_transpose(ri2h, ri2h);

    // Also store its inverse.
    mat3_invert(ri2h, rh2i);

    // Equatorial to ecliptic
    eraEcm06(DJM0, obs->tt, re2i);
    mat3_invert(re2i, ri2e);

    // For view matrix: use ecliptic in barycentric mode, else use ri2h
    if (obs->barycentric_mode) {
        // Barycentric mode - use ecliptic frame (ri2e) as "horizon"
        // Ecliptic plane becomes horizontal, ecliptic north pole becomes "up"
        mat3_mul(ro2v, ri2e, ri2v);  // ICRF → ecliptic → ro2v → view
    } else {
        mat3_mul(ro2v, ri2h, ri2v);  // ICRF → horizontal → ro2v → view
    }

    // ICRF to view (ignoring refraction).
    mat3_transpose(astrom->bpn, rc2v);
    if (obs->barycentric_mode) {
        // Barycentric mode - use ecliptic frame (ri2e) as "horizon"
        mat3_mul(ri2e, rc2v, rc2v);  // CIRS → ecliptic → view
    } else {
        mat3_mul(ri2h, rc2v, rc2v);  // CIRS → horizontal → view
    }
    mat3_mul(ro2v, rc2v, rc2v);  // Apply ro2v in both modes

    // Copy all
    mat3_copy(ro2v, obs->ro2v);
    mat3_invert(obs->ro2v, obs->rv2o);
    mat3_copy(ri2h, obs->ri2h);
    mat3_copy(rh2i, obs->rh2i);
    mat3_copy(ri2v, obs->ri2v);
    mat3_copy(ri2e, obs->ri2e);
    mat3_copy(re2i, obs->re2i);
    mat3_copy(rc2v, obs->rc2v);
}

static void observer_compute_hash(const observer_t *obs, uint64_t* hash_partial,
                                  uint64_t* hash)
{
    uint32_t v = 1;
    #define H(a) v = hash_xor(v, (const char*)&obs->a, sizeof(obs->a))
    H(elong);
    H(phi);
    H(hm);
    H(horizon);
    H(pressure);
    H(space);
    H(barycentric_mode);
    *hash_partial = v;
    H(ro2m);
    H(pitch);
    H(yaw);
    H(roll);
    H(view_offset_alt);
    H(tt);
    if (obs->space) H(obs_pvg);
    #undef H
    *hash = v;
}

static void correct_speed_of_light(double pv[2][3]) {
    double ldt = vec3_norm(pv[0]) * DAU2M / LIGHT_YEAR_IN_METER * DJY;
    vec3_addk(pv[0], pv[1], -ldt, pv[0]);
}


static void update_nutation_precession_mat(observer_t *obs)
{
    // XXX: we can maybe optimize this, since eraPn00a is very slow!
    double dpsi, deps, epsa, rb[3][3], rp[3][3], rbp[3][3], rn[3][3],
           rbpn[3][3];
    eraPn00a(DJM0, obs->tt, &dpsi, &deps, &epsa, rb, rp, rbp, rn, rbpn);
    mat3_mul(rn, rp, obs->rnp);

}

static void observer_update_fast(observer_t *obs)
{
    double dut1, theta, pvg[2][3];

    // Compute UT1 and UTC time.
    if (!obs->space && obs->last_update == obs->tt) goto end;

    obs->utc = tt2utc(obs->tt, &dut1);
    obs->ut1 = obs->utc + dut1 / ERFA_DAYSEC;

    if (!obs->space)
        eraAper13(DJM0, obs->ut1, &obs->astrom);
    eraPvu(obs->tt - obs->last_update, obs->earth_pvh, obs->earth_pvh);
    eraPvu(obs->tt - obs->last_update, obs->earth_pvb, obs->earth_pvb);

    if (!obs->space) {
        // Update observer geocentric position obs_pvg. We can't use eraPvu
        // here as the movement is a rotation about the earth center and can't
        // be approximated by a linear velocity  on a 24h time span
        theta = eraEra00(DJM0, obs->ut1);
        eraPvtob(obs->elong, obs->phi, obs->hm, 0, 0, 0, theta, obs->obs_pvg);
        // Rotate from CIRS to ICRF
        eraTrxp(obs->astrom.bpn, obs->obs_pvg[0], obs->obs_pvg[0]);
        eraTrxp(obs->astrom.bpn, obs->obs_pvg[1], obs->obs_pvg[1]);
        // Set pos back in AU
        eraSxp(DM2AU, obs->obs_pvg[0], obs->obs_pvg[0]);
        // Set speed back in AU / day
        eraSxp(ERFA_DAYSEC * DM2AU, obs->obs_pvg[1], obs->obs_pvg[1]);
    } else {
        vec3_mul(DAU2M, obs->obs_pvg[0], pvg[0]);
        vec3_mul(DAU2M / ERFA_DAYSEC, obs->obs_pvg[1], pvg[1]);
        eraApcs(DJM0, obs->tt, pvg, obs->earth_pvb, obs->earth_pvh[0],
                &obs->astrom);
    }

    // Update animation if active (uses core->clock for real-world time)
    if (obs->animating) {
        double t = (core->clock - obs->anim_start_clock) / obs->anim_duration;
        if (t >= 1.0) {
            // Animation complete
            eraCpv(obs->anim_target_pv, obs->barycentric_pv);
            obs->animating = false;
            LOG_E("Animation complete at clock=%.2f, final pos=[%.6f, %.6f, %.6f]",
                  core->clock, obs->barycentric_pv[0][0], obs->barycentric_pv[0][1],
                  obs->barycentric_pv[0][2]);
        } else {
            // Smooth interpolation (ease-in-out cubic)
            t = t < 0.5 ? 4 * t * t * t : 1 - pow(-2 * t + 2, 3) / 2;
            // Linear interpolation of position and velocity
            vec3_mix(obs->anim_start_pv[0], obs->anim_target_pv[0], t, obs->barycentric_pv[0]);
            vec3_mix(obs->anim_start_pv[1], obs->anim_target_pv[1], t, obs->barycentric_pv[1]);

            // Debug: log progress every 1%
            double total_dist_au = vec3_dist(obs->anim_start_pv[0], obs->anim_target_pv[0]);
            double current_percent = t * 100.0;

            // Log every 1% milestone
            if (floor(current_percent) > floor(obs->anim_last_logged_percent)) {
                double traveled[3];
                vec3_sub(obs->barycentric_pv[0], obs->anim_start_pv[0], traveled);
                double traveled_au = vec3_norm(traveled);

                // Format as light-years if distance > 1 LY, otherwise AU
                if (total_dist_au > 63241.0) {  // 1 LY ≈ 63241 AU
                    double traveled_ly = traveled_au / 63241.0;
                    double total_ly = total_dist_au / 63241.0;
                    LOG_E("Animation progress: %.2f ly / %.2f ly (%.0f%%), pos=[%.0f, %.0f, %.0f] AU",
                          traveled_ly, total_ly, current_percent,
                          obs->barycentric_pv[0][0], obs->barycentric_pv[0][1], obs->barycentric_pv[0][2]);
                } else {
                    LOG_E("Animation progress: %.0f AU / %.0f AU (%.0f%%), pos=[%.0f, %.0f, %.0f] AU",
                          traveled_au, total_dist_au, current_percent,
                          obs->barycentric_pv[0][0], obs->barycentric_pv[0][1], obs->barycentric_pv[0][2]);
                }
            }
            obs->anim_last_logged_percent = current_percent;
        }
        obs->barycentric_mode = true;
    }

    // Compute the observer's barycentric position
    if (obs->barycentric_mode) {
        // Use user-set barycentric position directly
        eraCpv(obs->barycentric_pv, obs->obs_pvb);
        // Calculate obs_pvg backwards for compatibility with frames.c/satellites.c
        eraPvmpv(obs->obs_pvb, obs->earth_pvb, obs->obs_pvg);
    } else {
        // Normal geocentric mode
        eraPvppv(obs->earth_pvb, obs->obs_pvg, obs->obs_pvb);
    }

end:
    update_matrices(obs);
    eraPvmpv(obs->earth_pvb, obs->earth_pvh, obs->sun_pvb);
    // Compute sun's apparent position in observer reference frame
    eraPvmpv(obs->sun_pvb, obs->obs_pvb, obs->sun_pvo);
    // Correct in one shot space motion, annual & diurnal abberrations
    correct_speed_of_light(obs->sun_pvo);
}

static void observer_update_full(observer_t *obs)
{
    double dut1, r[3][3], x, y, theta, s, sp, pvg[2][3];

    // Compute UT1 and UTC time.
    if (obs->last_update != obs->tt) {
        obs->utc = tt2utc(obs->tt, &dut1);
        obs->ut1 = obs->utc + dut1 / ERFA_DAYSEC;
    }

    // This is similar to a single call to eraApco13, except we handle
    // the time conversion ourself, since erfa doesn't support dates
    // before year -4800.
    eraPnm06a(DJM0, obs->tt, r); // equinox based BPN matrix.
    eraBpn2xy(r, &x, &y); // Extract CIP X,Y.
    s = eraS06(DJM0, obs->tt, x, y); // Obtain CIO locator s.
    // XXX: should be obs->ut1 here!  But it break the unit tests for now.
    theta = eraEra00(DJM0, obs->utc); // Earth rotation angle.
    sp = eraSp00(DJM0, obs->tt); // TIO locator s'.

    eraEpv00(DJM0, obs->tt, obs->earth_pvh, obs->earth_pvb);

    if (!obs->space) {
        eraApco(DJM0, obs->tt, obs->earth_pvb, obs->earth_pvh[0], x, y, s,
                theta, obs->elong, obs->phi, obs->hm, 0, 0, sp, 0, 0,
                &obs->astrom);
    } else {
        vec3_mul(DAU2M, obs->obs_pvg[0], pvg[0]);
        vec3_mul(DAU2M / ERFA_DAYSEC, obs->obs_pvg[1], pvg[1]);
        eraApcs(DJM0, obs->tt, pvg, obs->earth_pvb, obs->earth_pvh[0],
                &obs->astrom);
    }
    obs->eo = eraEors(r, s); // Equation of origins.

    // Update earth position.
    vec3_copy(obs->astrom.eb, obs->obs_pvb[0]);
    vec3_mul(ERFA_DC, obs->astrom.v, obs->obs_pvb[1]);
    if (!obs->space)
        eraPvmpv(obs->obs_pvb, obs->earth_pvb, obs->obs_pvg);
    // Update refraction constants.
    refraction_prepare(obs->pressure, 15, 0.5, &obs->refa, &obs->refb);
    update_nutation_precession_mat(obs);

    update_matrices(obs);
    eraPvmpv(obs->earth_pvb, obs->earth_pvh, obs->sun_pvb);
    // Compute sun's apparent position in observer reference frame
    eraPvmpv(obs->sun_pvb, obs->obs_pvb, obs->sun_pvo);
    // Correct in one shot space motion, annual & diurnal abberrations
    correct_speed_of_light(obs->sun_pvo);
}

EMSCRIPTEN_KEEPALIVE
void observer_update(observer_t *obs, bool fast)
{
    uint64_t hash, hash_partial;

    observer_compute_hash(obs, &hash_partial, &hash);
    // Check if we have computed accurate positions already
    if (hash == obs->hash)
        return;

    // Check if we have computed 'fast' positions already
    if (fast) {
        // Add one to the hash for the fast update hash value.
        hash++;
        if (hash == obs->hash)
            return;
        if (    hash_partial != obs->hash_partial ||
                fabs(obs->last_accurate_update - obs->tt) >= 1.001)
            fast = false;
    }

    if (fast)
        observer_update_fast(obs);
    else
        observer_update_full(obs);

    obs->last_update = obs->tt;
    obs->hash_partial = hash_partial;
    obs->hash = hash;
    if (!fast)
        obs->last_accurate_update = obs->tt;
}

static int observer_init(obj_t *obj, json_value *args)
{
    observer_t*  obs = (observer_t*)obj;
    mat3_set_identity(obs->ro2m);
    obs->tracking_target = NULL;
    obs->updating_tracking = false;
    observer_compute_hash(obs, &obs->hash_partial, &obs->hash);
    return 0;
}

static obj_t *observer_clone(const obj_t *obj)
{
    observer_t *ret;
    ret = (observer_t*)obj_create("observer", NULL);
    // Copy all except obj attributes.
    memcpy(((char*)ret) + sizeof(obj_t), ((const char*)obj) + sizeof(obj_t),
           sizeof(*ret) - sizeof(obj_t));
    return &ret->obj;
}

static void on_utc_changed(obj_t *obj, const attribute_t *attr)
{
    // Sync TT.
    observer_t *obs = (observer_t*)obj;
    obs->tt = utc2tt(obs->utc);
    module_changed(obj, "tt");
}

static void on_tt_changed(obj_t *obj, const attribute_t *attr)
{
    // Sync UTC.
    observer_t *obs = (observer_t*)obj;
    double dut1;
    obs->utc = tt2utc(obs->tt, &dut1);
    obs->ut1 = obs->utc + dut1 / ERFA_DAYSEC;
    module_changed(obj, "utc");
}

EMSCRIPTEN_KEEPALIVE
void observer_release_tracking(observer_t *obs)
{
    if (obs->tracking_target) {
        LOG_E("Releasing tracking target");
        obj_release(obs->tracking_target);
        obs->tracking_target = NULL;
    }
}

static void on_orientation_changed(obj_t *obj, const attribute_t *attr)
{
    // User manually changed orientation - stop tracking
    observer_t *obs = (observer_t*)obj;

    // Only stop tracking if we're actually tracking AND this isn't a tracking update
    if (obs->tracking_target && !obs->updating_tracking) {
        LOG_E("User manually changed %s - stopping tracking", attr->name);
        observer_release_tracking(obs);
    }
}

bool observer_is_uptodate(const observer_t *obs, bool fast)
{
    uint64_t hash, hash_partial;
    observer_compute_hash(obs, &hash_partial, &hash);
    if (hash == obs->hash) return true;
    if (fast && (hash + 1 == obs->hash)) return true;
    return false;
}

// Expose azalt vector to js.
static json_value *observer_get_azalt(obj_t *obj, const attribute_t *attr,
                                 const json_value *args)
{
    observer_t *obs = (observer_t*)obj;
    double v[3];
    vec3_from_sphe(obs->yaw, obs->pitch, v);
    return args_value_new(TYPE_V3, v);
}

// Expose/set obs_pvg position vector to/from js.
static json_value *observer_obs_pvg_pos_fn(obj_t *obj, const attribute_t *attr,
                                            const json_value *args)
{
    observer_t *obs = (observer_t*)obj;
    double v[3];
    // If args provided, this is a setter
    if (args && args->u.array.length) {
        args_get(args, TYPE_V3, v);
        vec3_copy(v, obs->obs_pvg[0]);
        return NULL;
    }
    // Otherwise, getter
    return args_value_new(TYPE_V3, obs->obs_pvg[0]);
}

// Expose/set obs_pvg velocity vector to/from js.
static json_value *observer_obs_pvg_vel_fn(obj_t *obj, const attribute_t *attr,
                                            const json_value *args)
{
    observer_t *obs = (observer_t*)obj;
    double v[3];
    // If args provided, this is a setter
    if (args && args->u.array.length) {
        args_get(args, TYPE_V3, v);
        vec3_copy(v, obs->obs_pvg[1]);
        return NULL;
    }
    // Otherwise, getter
    return args_value_new(TYPE_V3, obs->obs_pvg[1]);
}

// Expose/set barycentric position vector to/from js.
static json_value *observer_barycentric_position_fn(obj_t *obj, const attribute_t *attr,
                                                     const json_value *args)
{
    observer_t *obs = (observer_t*)obj;
    double v[3];
    // If args provided, this is a setter
    if (args && args->u.array.length) {
        args_get(args, TYPE_V3, v);
        // Check if setting to null/zero to disable barycentric mode
        if (vec3_norm2(v) < 1e-20) {
            obs->barycentric_mode = false;
            obs->space = false;  // Return to Earth surface mode
            // Re-enable atmosphere and landscape on Earth
            obj_t *atm = core_get_module("atmosphere");
            obj_t *ls = core_get_module("landscapes");
            if (atm) obj_set_attr(atm, "visible", true);
            if (ls) obj_set_attr(ls, "visible", true);
            // Restore mount frame to observed (horizon-based)
            core->mount_frame = FRAME_OBSERVED;
            core_update_mount(0);
        } else {
            vec3_copy(v, obs->barycentric_pv[0]);
            obs->barycentric_mode = true;
            obs->space = true;  // Enable space mode for correct aberration
            // Disable Earth-specific features in space
            obj_t *atm = core_get_module("atmosphere");
            obj_t *ls = core_get_module("landscapes");
            obj_t *lines = core_get_module("lines");
            if (atm) obj_set_attr(atm, "visible", false);
            if (ls) obj_set_attr(ls, "visible", false);
            // Disable horizon-based grids (azimuthal, meridian)
            if (lines) {
                obj_t *azimuthal = module_get_child(lines, "azimuthal");
                obj_t *meridian = module_get_child(lines, "meridian");
                if (azimuthal) obj_set_attr(azimuthal, "visible", false);
                if (meridian) obj_set_attr(meridian, "visible", false);
            }
            // Switch to ecliptic mount frame for proper roll/pitch/yaw
            core->mount_frame = FRAME_ECLIPTIC;
            core_update_mount(0);
        }
        return NULL;
    }
    // Otherwise, getter - return position or null if geocentric
    if (obs->barycentric_mode) {
        return args_value_new(TYPE_V3, obs->barycentric_pv[0]);
    }
    return json_null_new();
}

// Expose/set barycentric velocity vector to/from js.
static json_value *observer_barycentric_velocity_fn(obj_t *obj, const attribute_t *attr,
                                                     const json_value *args)
{
    observer_t *obs = (observer_t*)obj;
    double v[3];
    // If args provided, this is a setter
    if (args && args->u.array.length) {
        args_get(args, TYPE_V3, v);
        vec3_copy(v, obs->barycentric_pv[1]);
        obs->barycentric_mode = true;
        obs->space = true;  // Enable space mode for correct aberration
        // Disable Earth-specific features in space
        obj_t *atm = core_get_module("atmosphere");
        obj_t *ls = core_get_module("landscapes");
        obj_t *lines = core_get_module("lines");
        if (atm) obj_set_attr(atm, "visible", false);
        if (ls) obj_set_attr(ls, "visible", false);
        // Disable horizon-based grids (azimuthal, meridian)
        if (lines) {
            obj_t *azimuthal = module_get_child(lines, "azimuthal");
            obj_t *meridian = module_get_child(lines, "meridian");
            if (azimuthal) obj_set_attr(azimuthal, "visible", false);
            if (meridian) obj_set_attr(meridian, "visible", false);
        }
        // Switch to ecliptic mount frame for proper roll/pitch/yaw
        core->mount_frame = FRAME_ECLIPTIC;
        core_update_mount(0);
        return NULL;
    }
    // Otherwise, getter - return velocity or null if geocentric
    if (obs->barycentric_mode) {
        return args_value_new(TYPE_V3, obs->barycentric_pv[1]);
    }
    return json_null_new();
}

// Maximum safe travel distance to avoid numerical precision issues
// and prevent system hangs from extreme parallax/rendering calculations
#define MAX_SAFE_TRAVEL_LY 100.0   // Conservative safe limit (no precision loss)
#define MAX_HARD_TRAVEL_LY 1000.0  // Hard limit - beyond this causes hangs/precision loss

// C function to start animation - called from JavaScript wrapper
EMSCRIPTEN_KEEPALIVE
void observer_observe_from_object_ptr(observer_t *obs, obj_t *target, double duration_sec, double altitude_km, bool track)
{
    double pvo[2][4];

    if (!target || !obs) {
        LOG_E("observer_observe_from_object_ptr: null pointer");
        return;
    }

    // Get object's position/velocity
    // Try geometric position first (for planets), fall back to apparent (for stars)
    int ret = obj_get_info(target, obs, INFO_PVO_GEOMETRIC, pvo);
    if (ret != 0) {
        // Fall back to INFO_PVO for stars and other objects
        ret = obj_get_info(target, obs, INFO_PVO, pvo);
        if (ret != 0) {
            LOG_E("observeFromObject: failed to get object position");
            return;
        }
        LOG_E("observeFromObject: using apparent position (INFO_PVO) - target may be a star");
    }

    // Store Earth-centric distance before converting to barycentric (needed for radius calculation)
    double earth_centric_distance = vec3_norm(pvo[0]);

    // Check if object is at infinity (stars, distant objects)
    if (pvo[0][3] == 0.0) {
        // Object is at infinity - pvo[0] is a unit direction vector
        // Need actual distance to position observer
        double distance;
        ret = obj_get_info(target, obs, INFO_DISTANCE, &distance);
        if (ret != 0 || distance <= 0 || isnan(distance)) {
            LOG_E("observeFromObject: cannot observe from object at infinity without valid distance (ret=%d, dist=%.2f)",
                  ret, distance);
            return;
        }

        double distance_ly = distance / 63241.077;
        LOG_E("observeFromObject: object at infinity, distance = %.2f AU (%.2f light-years)",
              distance, distance_ly);

        // Check distance limits to avoid numerical precision issues
        if (distance_ly > MAX_HARD_TRAVEL_LY) {
            LOG_E("ERROR: Distance %.0f LY exceeds hard limit of %.0f LY - precision loss will occur",
                  distance_ly, MAX_HARD_TRAVEL_LY);
            LOG_E("Cannot travel to objects beyond %.0f light-years due to floating point precision limits",
                  MAX_HARD_TRAVEL_LY);
            return;
        } else if (distance_ly > MAX_SAFE_TRAVEL_LY) {
            LOG_E("WARNING: Distance %.0f LY exceeds recommended limit of %.0f LY - may cause rendering issues",
                  distance_ly, MAX_SAFE_TRAVEL_LY);
            LOG_E("Proceeding anyway, but expect potential visual artifacts");
        }

        // Scale unit direction by distance to get actual position
        vec3_mul(distance, pvo[0], pvo[0]);
        pvo[0][3] = 1.0;  // No longer at infinity
        earth_centric_distance = distance;
        // Velocity stays zero for stars
    }

    // Convert to barycentric by adding Earth's barycentric position
    // pvo is currently Earth-centric, convert to barycentric
    vec3_add(pvo[0], obs->earth_pvb[0], pvo[0]);
    vec3_add(pvo[1], obs->earth_pvb[1], pvo[1]);

    // Validate resulting position
    double pos_magnitude = vec3_norm(pvo[0]);
    if (isnan(pos_magnitude)) {
        LOG_E("observeFromObject: invalid position after barycentric conversion (magnitude=%f)", pos_magnitude);
        return;
    }

    // Final sanity check on total distance (catches edge cases)
    double final_distance_ly = pos_magnitude / 63241.077;
    if (final_distance_ly > MAX_HARD_TRAVEL_LY) {
        LOG_E("ERROR: Final barycentric position %.0f LY exceeds hard limit", final_distance_ly);
        return;
    }

    LOG_E("observeFromObject: target barycentric pos = [%.6f, %.6f, %.6f] AU (magnitude=%.0f AU, %.2f LY)",
          pvo[0][0], pvo[0][1], pvo[0][2], pos_magnitude, final_distance_ly);

    // Position observer at specified altitude above the surface
    // For objects without radius info (stars), position at their center
    double radius_au = 0;
    int has_radius = 0;

    // Try INFO_PHYSICAL_RADIUS first (planets)
    if (obj_get_info(target, obs, INFO_PHYSICAL_RADIUS, &radius_au) == 0 && radius_au > 0) {
        has_radius = 1;
    } else {
        // Fall back to INFO_RADIUS (angular) for asteroids
        double angular_radius;
        if (obj_get_info(target, obs, INFO_RADIUS, &angular_radius) == 0 && angular_radius > 0) {
            // Convert from angular radius to physical radius
            // angular_radius is in radians, earth_centric_distance is in AU
            if (earth_centric_distance > 0) {
                radius_au = angular_radius * earth_centric_distance;
                has_radius = 1;
                LOG_E("observeFromObject: using INFO_RADIUS, angular=%.10f rad, distance=%.6f AU",
                      angular_radius, earth_centric_distance);
            }
        }
    }

    if (has_radius) {
        LOG_E("observeFromObject: target radius = %.10f AU (%.1f km)",
              radius_au, radius_au * 149597870.7);
        LOG_E("observeFromObject: sun barycentric = [%.6f, %.6f, %.6f]",
              obs->sun_pvb[0][0], obs->sun_pvb[0][1], obs->sun_pvb[0][2]);

        // Convert altitude to AU: km / 149,597,870.7 km/AU
        // Default to 10,000 km if not specified
        if (altitude_km <= 0) altitude_km = 10000.0;
        const double altitude_au = altitude_km / 149597870.7;
        double surface_distance_au = radius_au + altitude_au;

        // Direction: from planet center towards the Sun (to view day side)
        double direction[3];
        vec3_sub(obs->sun_pvb[0], pvo[0], direction);

        double dir_magnitude = vec3_norm(direction);
        LOG_E("observeFromObject: direction magnitude before norm = %.6f AU", dir_magnitude);

        // If sun direction is zero (observing the Sun itself), use ecliptic north
        if (vec3_norm2(direction) < 1e-20) {
            direction[0] = 0;
            direction[1] = 0;
            direction[2] = 1;
        }

        vec3_normalize(direction, direction);
        LOG_E("observeFromObject: normalized direction = [%.6f, %.6f, %.6f]",
              direction[0], direction[1], direction[2]);

        // Calculate offset from planet center
        double offset[3];
        vec3_mul(surface_distance_au, direction, offset);

        LOG_E("observeFromObject: offset = [%.10f, %.10f, %.10f] AU (magnitude %.10f AU)",
              offset[0], offset[1], offset[2], vec3_norm(offset));

        // Add offset to target position
        vec3_add(pvo[0], offset, pvo[0]);
        // Velocity remains the same (moving with the planet)

        LOG_E("observeFromObject: final position = [%.6f, %.6f, %.6f] AU",
              pvo[0][0], pvo[0][1], pvo[0][2]);
    } else {
        LOG_E("observeFromObject: no physical radius available, positioning at object center");
    }

    LOG_E("observeFromObject: duration = %.2f sec, clock = %.2f",
          duration_sec, core->clock);

    // If duration is 0, jump instantly without animation
    if (duration_sec <= 0) {
        LOG_E("observeFromObject: instant jump (no animation)");
        vec3_copy(pvo[0], obs->barycentric_pv[0]);
        vec3_copy(pvo[1], obs->barycentric_pv[1]);
        obs->animating = false;
    } else {
        // Set up animation (duration is in real-world seconds)
        obs->anim_start_clock = core->clock;
        obs->anim_duration = duration_sec;

        // Current position as start (or Earth if not in barycentric mode)
        if (obs->barycentric_mode) {
            eraCpv(obs->barycentric_pv, obs->anim_start_pv);
            LOG_D("observeFromObject: starting from barycentric [%.6f, %.6f, %.6f]",
                  obs->anim_start_pv[0][0], obs->anim_start_pv[0][1], obs->anim_start_pv[0][2]);
        } else {
            // Start from Earth's position
            vec3_set(obs->anim_start_pv[0], 0, 0, 0);
            vec3_set(obs->anim_start_pv[1], 0, 0, 0);
            LOG_D("observeFromObject: starting from Earth [0, 0, 0]");
        }

        // Target is the object's position
        vec3_copy(pvo[0], obs->anim_target_pv[0]);
        vec3_copy(pvo[1], obs->anim_target_pv[1]);

        // Start animation
        obs->animating = true;
        obs->anim_last_logged_percent = 0.0;  // Reset progress logging
    }

    obs->barycentric_mode = true;
    obs->space = true;  // Enable space mode for correct aberration

    // Reset roll to align with ecliptic plane
    obs->roll = 0;
    obs->pitch = 0;
    obs->yaw = 0;

    // Disable Earth-specific features in space
    obj_t *atm = core_get_module("atmosphere");
    obj_t *ls = core_get_module("landscapes");
    obj_t *lines = core_get_module("lines");
    if (atm) obj_set_attr(atm, "visible", false);
    if (ls) obj_set_attr(ls, "visible", false);
    // Disable horizon-based grids (azimuthal, meridian)
    if (lines) {
        obj_t *azimuthal = module_get_child(lines, "azimuthal");
        obj_t *meridian = module_get_child(lines, "meridian");
        if (azimuthal) obj_set_attr(azimuthal, "visible", false);
        if (meridian) obj_set_attr(meridian, "visible", false);
    }
    // Switch to ecliptic mount frame for proper roll/pitch/yaw
    core->mount_frame = FRAME_ECLIPTIC;
    // Force immediate mount update so ro2m gets updated before next matrix computation
    core_update_mount(0);

    // Set up tracking if requested
    if (track) {
        obj_retain(target);  // Retain target so it doesn't get freed
        observer_release_tracking(obs);
        obs->tracking_target = target;
        LOG_E("observeFromObject: tracking enabled for target");
    } else {
        observer_release_tracking(obs);
        LOG_E("observeFromObject: tracking disabled");
    }

    LOG_D("observeFromObject: animation started, orientation reset to ecliptic!");
}


static obj_klass_t observer_klass = {
    .id = "observer",
    .size = sizeof(observer_t),
    .flags = OBJ_IN_JSON_TREE,
    .init = observer_init,
    .clone = observer_clone,
    .attributes = (attribute_t[]) {
        PROPERTY(longitude, TYPE_ANGLE, MEMBER(observer_t, elong)),
        PROPERTY(latitude, TYPE_ANGLE, MEMBER(observer_t, phi)),
        PROPERTY(elevation, TYPE_FLOAT, MEMBER(observer_t, hm)),
        PROPERTY(tt, TYPE_MJD, MEMBER(observer_t, tt),
                 .on_changed = on_tt_changed),
        PROPERTY(utc, TYPE_MJD, MEMBER(observer_t, utc),
                 .on_changed = on_utc_changed),
        PROPERTY(pitch, TYPE_ANGLE, MEMBER(observer_t, pitch),
                 .on_changed = on_orientation_changed),
        PROPERTY(yaw, TYPE_ANGLE, MEMBER(observer_t, yaw),
                 .on_changed = on_orientation_changed),
        PROPERTY(roll, TYPE_ANGLE, MEMBER(observer_t, roll),
                 .on_changed = on_orientation_changed),
        PROPERTY(view_offset_alt, TYPE_ANGLE,
                 MEMBER(observer_t, view_offset_alt)),
        PROPERTY(azalt, TYPE_V3, .fn = observer_get_azalt),
        PROPERTY(space, TYPE_BOOL, MEMBER(observer_t, space)),
        PROPERTY(obs_pvg_pos, TYPE_V3, .fn = observer_obs_pvg_pos_fn),
        PROPERTY(obs_pvg_vel, TYPE_V3, .fn = observer_obs_pvg_vel_fn),
        PROPERTY(barycentricPosition, TYPE_V3, .fn = observer_barycentric_position_fn),
        PROPERTY(barycentricVelocity, TYPE_V3, .fn = observer_barycentric_velocity_fn),
        {}
    },
};
OBJ_REGISTER(observer_klass)

