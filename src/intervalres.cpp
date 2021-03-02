#include "intervalres.h"
#include "helper.h"

IntervalRes::IntervalRes()
{
}


IntervalRes::IntervalRes(std::string _id, std::string _chr, int _start, int _end) {
    id = _id;
    chr = _chr;
    start = _start;
    end = _end;
}

IntervalRes::IntervalRes(std::string _desc) {
    set_data(_desc);
}

IntervalRes::IntervalRes(std::string _id, std::string _desc) {
    id = _id;
    set_data(_desc);
}


void IntervalRes::set_data(std::string _desc) {

    // ID:CHR:START..END  or CHR:START..END
    std::vector<std::string> fields = tokenize_string(_desc, ':');

	
    if(fields.size()==3) {
        id = fields.at(0);
    	chr = fields.at(1);
    }
    else if (fields.size()!=2) {
        // Can happen!
        return;
    }
    else {
    	chr = fields.at(0);
    }

    std::vector<std::string> loc = tokenize_string( fields.at(fields.size()-1), '.');
    if(loc.size()==2) {
        str2int(loc.at(0), start);
        str2int(loc.at(1), end);
    }
    else if(loc.size()==1) {
        str2int(loc.at(0), start);
        end = -1;
    }
}

std::string IntervalRes::get_loc() {
    if(end == -1) {
        // SNP
        return chr + ":" + int2str(start);
    }
    else {
        return chr + ":" + int2str(start) + ".." + int2str(end);
    }
}
