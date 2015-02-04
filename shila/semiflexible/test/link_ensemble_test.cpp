#define BOOST_TEST_MODULE link_ensemble_test

#include "Link.h"
#include "actin_ensemble.h"
#include "link_ensemble.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    actin * a1 = new actin(-1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin( 1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a1, 0);
    ae.add_monomer(a2, 0);
   
    link_ensemble l = link_ensemble();
} 

BOOST_AUTO_TEST_CASE( link_walk_test )
{
    double lnk_len = 1;
    double kl = 100;
    double kb = 2;
    std::string col = "yellow";
    
    actin * a1 = new actin(-1, 2, 0, 1, 0, 0, 0, 0, 0);
    actin * a2 = new actin( 1, 2, 0, 1, 0, 0, 0, 0, 0);
    
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a1, 0);
    ae.add_monomer(a2, 0);
    link_ensemble l = link_ensemble();
    
    ae.connect_polymers(&l, lnk_len, kl, kb, col);
    l.link_walk();

}
BOOST_AUTO_TEST_CASE( friction_test )
{
}

BOOST_AUTO_TEST_CASE( force_test)
{
}

BOOST_AUTO_TEST_CASE( direction_test)
{
}

BOOST_AUTO_TEST_CASE( start_end_test)
{
}

BOOST_AUTO_TEST_CASE( get_intpoint_test)
{
}
// EOF
