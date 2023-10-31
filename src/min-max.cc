#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{
    double a = (p2 - p3).two_norm();
    double b = (p3 - p1).two_norm();
    double c = (p1 - p2).two_norm();

    double cos_beta = (a * a + c * c - b * b) / (2 * a * c);
    double rez = std::acos(std::min(1.0, std::max(-1.0, cos_beta))) * 180.0 / M_PI;

    return rez;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using LeafGridView = GridType::LeafGridView;

    std::unique_ptr<GridType> nas_grid = Dune::GmshReader<GridType>::read(argv[1]);
    auto grid = nas_grid->leafGridView();

    double max_kut = 0;
    double min_kut = 180;
    int broj_elemenata = 0;

    for (auto const & element : elements(grid))
    {
        auto geom = element.geometry();
        auto p1 = geom.corner(0);
        auto p2 = geom.corner(1);
        auto p3 = geom.corner(2);

        double kut1 = kut(p1, p2, p3);
        double kut2 = kut(p2, p3, p1);
        double kut3 = kut(p3, p1, p2);

        max_kut = std::max({max_kut, kut1, kut2, kut3});
        min_kut = std::min({min_kut, kut1, kut2, kut3});
        broj_elemenata++;
    }

    std::cout << std::endl;
    std::cout << "Broj elemenata je : " << broj_elemenata << "." << std::endl;
    std::cout << "Minimalni kut je : " << min_kut << "." << std::endl;
    std::cout << "Maksimalni kut je : " << max_kut << "." << std::endl;

    Dune::VTKWriter<LeafGridView> vtkwriter(grid);
    vtkwriter.write("poluvijenac");

    return 0;
}
