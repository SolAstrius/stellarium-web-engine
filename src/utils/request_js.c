/* Stellarium Web Engine - Copyright (c) 2022 - Stellarium Labs SRL
 *
 * This program is licensed under the terms of the GNU AGPL v3, or
 * alternatively under a commercial licence.
 *
 * The terms of the AGPL v3 license can be found in the main directory of this
 * repository.
 */

#include "swe.h"

#ifdef __EMSCRIPTEN__

#include <emscripten/fetch.h>

#define MAX_NB  16 // Max number of concurrent requests.

struct request
{
    char        *url;
    emscripten_fetch_t *fetch;
    int         status_code;
    bool        done;
    void        *data;
    int         size;
};


static struct {
    int nb;     // Number of current running requests.
} g = {};

static bool url_has_extension(const char *str, const char *ext);

void request_init(const char *cache_dir)
{
    // Ignore the cache dir with emscripten.

    // Some checks that 'url_has_extension' works well enough.
    assert(url_has_extension("https://xyz.test.jpg", ".jpg"));
    assert(url_has_extension("http://xyz.test.jpg?xyz", ".jpg"));
    assert(url_has_extension("http://xyz.test.jpg#xyz", ".jpg"));
}

request_t *request_create(const char *url)
{
    request_t *req = calloc(1, sizeof(*req));
    req->url = strdup(url);
    return req;
}

int request_is_finished(const request_t *req)
{
    return req->done;
}

void request_delete(request_t *req)
{
    if (!req) return;
    if (req->fetch) {
        emscripten_fetch_close(req->fetch);
        g.nb--;
    }
    free(req->url);
    free(req->data);
    free(req);
}

// Very basic algo to check the extension of a url file.
// XXX: should probably use a regex.
static bool url_has_extension(const char *str, const char *ext)
{
    const char *s;
    char c;

    s = strcasestr(str, ext);
    if (!s) return false;
    c = s[strlen(ext)];
    if (c != '\0' && c != '?' && c != '#') return false;
    return true;
}

static bool could_be_str(const request_t *req)
{
    return !url_has_extension(req->url, ".jpeg") &&
           !url_has_extension(req->url, ".jpg") &&
           !url_has_extension(req->url, ".png") &&
           !url_has_extension(req->url, ".webp") &&
           !url_has_extension(req->url, ".eph");
}

static void downloadSucceeded(emscripten_fetch_t *fetch)
{
    request_t *req = (request_t*)fetch->userData;

    req->status_code = fetch->status;
    req->size = fetch->numBytes;

    // Copy the data from fetch buffer
    req->data = malloc(req->size + 1);
    memcpy(req->data, fetch->data, req->size);

    // Add null termination for text files
    if (could_be_str(req)) {
        ((char*)req->data)[req->size] = '\0';
    }

    req->done = true;
    req->fetch = NULL;
    emscripten_fetch_close(fetch);
    g.nb--;
}

static void downloadFailed(emscripten_fetch_t *fetch)
{
    request_t *req = (request_t*)fetch->userData;

    req->status_code = fetch->status ?: 499;
    req->done = true;
    req->fetch = NULL;
    emscripten_fetch_close(fetch);
    g.nb--;
}

const void *request_get_data(request_t *req, int *size, int *status_code)
{
    if (!req->done && !req->fetch && g.nb < MAX_NB) {
        emscripten_fetch_attr_t attr;
        emscripten_fetch_attr_init(&attr);
        strcpy(attr.requestMethod, "GET");
        attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
        attr.onsuccess = downloadSucceeded;
        attr.onerror = downloadFailed;
        attr.userData = req;

        req->fetch = emscripten_fetch(&attr, req->url);
        g.nb++;
    }
    if (size) *size = req->size;
    if (status_code) *status_code = req->status_code;
    return req->data;
}

void request_make_fresh(request_t *req)
{
}

#endif
