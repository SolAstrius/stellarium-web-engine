# Phase 1: Critical Fixes & Compatibility

## Goal
Fix breaking bugs from Emscripten 4.x while maintaining existing API.

## Changes Required

### 1.1 Fix Missing Export (obj.js:382-384)

**Current (BROKEN):**
```javascript
function stringToC(str) {
    let size = Module.lengthBytesUTF8(str);
    let ptr = Module._malloc(size + 1);
    Module.writeAsciiToMemory(str, ptr, false);  // ❌ Not exported!
    return ptr;
}
```

**Fixed:**
```javascript
function stringToC(str) {
    let size = Module.lengthBytesUTF8(str);
    let ptr = Module._malloc(size + 1);
    Module.stringToUTF8(str, ptr, size + 1);  // ✅ Works in Em 4.x
    return ptr;
}
```

### 1.2 Update SConstruct Exports

**Current:**
```python
extra_exported = [
    'GL',
    'UTF8ToString',
    'addFunction',      # Still needed for Phase 1
    'ccall',
    'cwrap',
    'getValue',
    'intArrayFromString',
    'lengthBytesUTF8',
    'removeFunction',
    'setValue',
    'stringToUTF8',
]
```

**Add (needed but missing):**
```python
extra_exported = [
    'GL',
    'UTF8ToString',
    'stringToUTF8',      # ← ADD THIS
    'addFunction',       # Keep for now
    'removeFunction',    # Keep for now
    'ccall',
    'cwrap',
    'getValue',
    'setValue',
    'intArrayFromString',
    'lengthBytesUTF8',
]
```

### 1.3 Verify All Wrapped Functions Have KEEPALIVE

Run this audit script:

```bash
# Extract all cwrapped functions
grep -h "Module.cwrap" src/js/*.js | \
  grep -oP "'[a-z_]+'" | sort -u > /tmp/wrapped_funcs.txt

# Check each has EMSCRIPTEN_KEEPALIVE
while read func; do
  func_clean=$(echo $func | tr -d "'")
  if ! grep -r "EMSCRIPTEN_KEEPALIVE" src/ -A1 | grep -q "$func_clean"; then
    echo "❌ MISSING KEEPALIVE: $func_clean"
  fi
done < /tmp/wrapped_funcs.txt
```

### 1.4 Fix Property Getters Returning Null

**Current issue in observer.c:389-411:**

The barycentric property functions return `json_null_new()` which becomes `null` in JS, causing:
```
TypeError: ret is null
```

**Fix:** Check for null in JavaScript before parsing:

```javascript
// In obj.js _call() method (around line 320):
SweObj.prototype._call = function(attr, arg) {
    if (arg === undefined || arg === null)
        arg = 0
    else
        arg = JSON.stringify(arg)
    var cret = obj_call_json_str(this.v, attr, arg)

    // ✅ ADD NULL CHECK:
    if (cret === 0) return null;

    var ret = Module.UTF8ToString(cret)
    Module._free(cret)
    if (!ret) return null;
    ret = JSON.parse(ret)
    if (!ret.swe_) return ret;
    if (ret.type === 'obj') {
        let v = parseInt(ret.v)
        return v ? new SweObj(v) : null;
    }
    return ret.v;
}
```

## Estimated Time
- **2-4 hours** of focused work
- **Low risk** - maintains compatibility

## Testing Checklist

```javascript
// After rebuild, test in console:

// 1. Test string conversion works
stel.getObj('NAME Mars')  // Should return object, not error

// 2. Test barycentric properties
stel.observer.barycentricPosition  // Should return null (not crash)

// 3. Test observeFromObject
const mars = stel.getObj('NAME Mars')
stel.observer.observeFromObject(mars, 5.0)  // Should work

// 4. Test other properties still work
stel.observer.longitude
stel.observer.latitude
stel.core.observer.tt
```

## Success Criteria
- ✅ No console errors on page load
- ✅ `stel.getObj()` works
- ✅ Properties return values without crashing
- ✅ `observeFromObject()` animates camera
