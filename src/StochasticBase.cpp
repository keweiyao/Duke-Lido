#include "StochasticBase.h"
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

template<size_t N>
StochasticBase<N>::StochasticBase(std::string Name, 
								boost::property_tree::ptree config):
_Name(Name)
{
	std::vector<std::string> strs, slots;
	boost::split(strs, _Name, boost::is_any_of("/"));
	std::string process_name = strs[0];
	std::string quantity_name = strs[1];
	std::cout << __func__ << " " << _Name << std::endl;
	auto tree = config.get_child(process_name+"."+quantity_name);
	std::string allslots = tree.get<std::string>("<xmlattr>.slots");
	boost::split(slots, allslots, boost::is_any_of(",") );

	std::vector<size_t> shape;
	std::vector<double> low, high;
	for(auto & v : slots){
		shape.push_back(tree.get<size_t>("N"+v));
		low.push_back(tree.get<double>("L"+v));
		high.push_back(tree.get<double>("H"+v));
	}
	
	StochasticBase<N>::_ZeroMoment = 
		std::make_shared<TableBase<scalar, N>>(Name+"/scalar", shape, low, high);
	StochasticBase<N>::_FirstMoment = 
		std::make_shared<TableBase<fourvec, N>>(Name+"/vector", shape, low, high);
	StochasticBase<N>::_SecondMoment = 
		std::make_shared<TableBase<tensor, N>>(Name+"/tensor", shape, low, high);
_ZeroMoment->SetTableValue({0, 0}, {10});
	_ZeroMoment->Save("table.h5");
	_FirstMoment->Save("table.h5");
	_SecondMoment->Save("table.h5");
}

template class StochasticBase<2>;

