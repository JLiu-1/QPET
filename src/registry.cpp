#include <stdexcept>

#include "qpet/backend.hpp"

namespace qpet {

// Each backend lives in its own translation unit, because the compressors' headers do not all coexist.
const Backend &sz3_backend();
const Backend &szo_backend();
const Backend &zfp_backend();
const Backend &sperr_backend();

namespace {

std::vector<const Backend *> backends() {
    return {
#ifdef QPET_WITH_SZ3
        &sz3_backend(),
#endif
#ifdef QPET_WITH_SZO
        &szo_backend(),
#endif
#ifdef QPET_WITH_ZFP
        &zfp_backend(),
#endif
#ifdef QPET_WITH_SPERR
        &sperr_backend(),
#endif
    };
}

}  // namespace

const Backend &get_backend(const std::string &name) {
    for (const Backend *b : backends())
        if (b->name() == name) return *b;
    throw std::invalid_argument("QPET: this build has no backend named \"" + name + "\"");
}

std::vector<std::string> available_backends() {
    std::vector<std::string> names;
    for (const Backend *b : backends()) names.push_back(b->name());
    return names;
}

}  // namespace qpet
