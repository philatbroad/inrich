#include "interval.h"
#include "helper.h"
#include "loaders.h"

std::string Interval::get_desc() {
    return name + ":chr"  + get_hum_chrcode(chr) + ":" + int2str(bp1) + ".." + int2str(bp2);
}
