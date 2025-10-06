/* Stellarium Web Engine - Copyright (c) 2022 - Stellarium Labs SRL
 *
 * This program is licensed under the terms of the GNU AGPL v3, or
 * alternatively under a commercial licence.
 *
 * The terms of the AGPL v3 license can be found in the main directory of this
 * repository.
 */

#include "swe.h"

#define D 0.7071067811865476 // Sqrt(2)/2

static struct {
    const char *text;          // Earth cardinal direction
    const char *text_ecliptic; // Ecliptic/space cardinal direction
    double pos[3];
} POINTS[] = {
    {"N" , "ECL N", { 1,  0,  0}},  // North → Ecliptic North
    {"E" , "ECL E", { 0,  1,  0}},  // East → Ecliptic East
    {"S" , "ECL S", {-1,  0,  0}},  // South → Ecliptic South
    {"W" , "ECL W", { 0, -1,  0}},  // West → Ecliptic West
    {"" , "VE", { 0,  0,  0}},      // Vernal Equinox (space only)
    {"" , "AE", { 0,  0,  0}},      // Autumnal Equinox (space only)

    {"NE", "ECL NE", { D,  D,  0}},
    {"SE", "ECL SE", {-D,  D,  0}},
    {"SW", "ECL SW", {-D, -D,  0}},
    {"NW", "ECL NW", { D, -D,  0}},
};


/*
 * Type: cardinal_t
 * Cardinal module.
 */
typedef struct cardinal {
    obj_t obj;
    fader_t         visible;
} cardinal_t;

static int cardinal_init(obj_t *obj, json_value *args)
{
    cardinal_t *c = (void*)obj;
    fader_init(&c->visible, true);
    return 0;
}

static int cardinal_update(obj_t *obj, double dt)
{
    cardinal_t *c = (void*)obj;
    return fader_update(&c->visible, dt);
}

static int cardinal_render(obj_t *obj, const painter_t *painter)
{
    int i;
    double size = 24;
    const cardinal_t *c = (const cardinal_t*)obj;
    double color[4] = {0.8, 0.4, 0.4, 0.8 * c->visible.value};
    double p[4], ecliptic_pos[3];
    const char *label;
    bool in_space = painter->obs->barycentric_mode;

    // Ecliptic cardinal directions in ecliptic frame, then converted to ICRF
    double ecl_dirs[6][3] = {
        {0, 0, 1},   // ECL N = ecliptic north pole (Z-axis in ecliptic)
        {0, 1, 0},   // ECL E = 90° along ecliptic (Y-axis in ecliptic)
        {0, 0, -1},  // ECL S = ecliptic south pole (-Z in ecliptic)
        {0, -1, 0},  // ECL W = 270° along ecliptic (-Y in ecliptic)
        {1, 0, 0},   // VE = vernal equinox, 0° ecliptic longitude (X-axis)
        {-1, 0, 0},  // AE = autumnal equinox, 180° ecliptic longitude (-X-axis)
    };

    if (c->visible.value <= 0) return 0;
    int num_points = in_space ? 6 : 4;  // Show VE/AE only in space
    for (i = 0; i < num_points; i++) {
        // In space, use ecliptic coordinates
        if (in_space) {
            // Convert ecliptic direction to ICRF
            mat3_mul_vec3(painter->obs->re2i, ecl_dirs[i], ecliptic_pos);

            if (painter_is_point_clipped_fast(painter, FRAME_ICRF,
                    ecliptic_pos, true))
                continue;
            convert_frame(painter->obs, FRAME_ICRF, FRAME_VIEW, true,
                          ecliptic_pos, p);
        } else {
            // On Earth, use horizon frame as normal
            if (painter_is_point_clipped_fast(painter, FRAME_OBSERVED,
                    POINTS[i].pos, true))
                continue;
            convert_frame(painter->obs, FRAME_OBSERVED, FRAME_VIEW, true,
                          POINTS[i].pos, p);
        }

        project_to_win(painter->proj, p, p);
        painter_t _painter = *painter;
        _painter.color[0] = color[0];
        _painter.color[1] = color[1];
        _painter.color[2] = color[2];
        _painter.color[3] = c->visible.value;
        _painter.lines.width = 4;
        paint_2d_ellipse(&_painter, NULL, 0, p, VEC(1, 1), NULL);

        // Use ecliptic labels in space, geographic labels on Earth
        label = in_space ? POINTS[i].text_ecliptic : POINTS[i].text;

        if (in_space) {
            labels_add_3d(sys_translate("gui", label), FRAME_ICRF,
                          ecliptic_pos, true, 0, size, color, 0,
                          ALIGN_CENTER | ALIGN_TOP, TEXT_BOLD, 0, NULL);
        } else {
            labels_add_3d(sys_translate("gui", label), FRAME_OBSERVED,
                          POINTS[i].pos, true, 0, size, color, 0,
                          ALIGN_CENTER | ALIGN_TOP, TEXT_BOLD, 0, NULL);
        }
    }
    return 0;
}


/*
 * Meta class declaration.
 */

static obj_klass_t cardinal_klass = {
    .id = "cardinals",
    .size = sizeof(cardinal_t),
    .flags = OBJ_IN_JSON_TREE | OBJ_MODULE,
    .render = cardinal_render,
    .init   = cardinal_init,
    .update = cardinal_update,
    .render_order = 50,
    .attributes = (attribute_t[]) {
        PROPERTY(visible, TYPE_BOOL, MEMBER(cardinal_t, visible.target)),
        {}
    },
};
OBJ_REGISTER(cardinal_klass)
