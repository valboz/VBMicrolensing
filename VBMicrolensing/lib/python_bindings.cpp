#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "VBMicrolensingLibrary.h"
#include <string>
#include <pybind11/functional.h>




namespace py = pybind11;

// Declaration of an instance to VBMicrolensing class. 
VBMicrolensing VBM;

PYBIND11_MODULE(VBMicrolensing, m) {
    py::options options;
    options.disable_function_signatures();

    py::class_<VBMicrolensing> vbm(m, "VBMicrolensing");
    vbm.def(py::init());
    // Settings

    vbm.def_readwrite("Tol", &VBMicrolensing::Tol,
        "Absolute accuracy goal.");
    vbm.def_readwrite("RelTol", &VBMicrolensing::RelTol,
        "Relative precision goal.");
    vbm.def_readwrite("a1", &VBMicrolensing::a1,
        "Linear limb darkening coefficient. I(r)=I(0)(1-a1(1-\\sqrt{1-r^2/\\rho^2}))");
    vbm.def_readwrite("a2", &VBMicrolensing::a2,
        "Second limb darkening coefficient.");
    vbm.def_readwrite("minannuli", &VBMicrolensing::minannuli,
        "Minimum number of annuli to calculate for limb darkening.");
    vbm.def_readwrite("NPcrit", &VBMicrolensing::NPcrit,
        "Number of points in critical curves.");
    vbm.def_readwrite("parallaxsystem", &VBMicrolensing::parallaxsystem,
        "0 for parallel-perpendicular, 1 for North-Eeast.");
    vbm.def_readwrite("t0_par_fixed", &VBMicrolensing::t0_par_fixed,
        "Set to 1 if you want to specify a constant t_{0,par}.");
    vbm.def_readwrite("t0_par", &VBMicrolensing::t0_par,
        "Reference time for parallax t_{0,par}. Only used if t0_par_fixed=1.");
    vbm.def_readwrite("satellite", &VBMicrolensing::satellite,
        "Specifies the satellite number for the next calculation \
                (0 for observations from the ground);.");
    vbm.def_readwrite("astrometry", &VBMicrolensing::astrometry,
        "Unlock astrometry centroid calculation.");
    vbm.def_readwrite("astrox1", &VBMicrolensing::astrox1,
        "The x component of the light centroid.");
    vbm.def_readwrite("astrox2", &VBMicrolensing::astrox2,
        "The y component of the light centroid.");
    vbm.def_readwrite("mass_luminosity_exponent", &VBMicrolensing::mass_luminosity_exponent,
        "Exponent for the mass-luminosity relation: L = M^q; default value is q=4.0");
    vbm.def_readwrite("mass_radius_exponent", &VBMicrolensing::mass_radius_exponent,
        "Exponent for the mass-radius relation: R = M^q; default value is q=0.89");
    vbm.def_readwrite("lens_mass_luminosity_exponent", &VBMicrolensing::lens_mass_luminosity_exponent,
        "Exponent for the mass-luminosity relation for the lens: L = M^q; default value is q=4.0");
    vbm.def_readwrite("turn_off_secondary_source", &VBMicrolensing::turn_off_secondary_source,
        "Flux of secondary source is set to zero.");
    vbm.def_readwrite("turn_off_secondary_lens", &VBMicrolensing::turn_off_secondary_lens,
        "Flux of secondary lens is set to zero.");
    vbm.def_readwrite("corrquad", &VBMicrolensing::corrquad,
        "Quadrupole test.");
    vbm.def_readwrite("corrquad2", &VBMicrolensing::corrquad2,
        "Ghost image test.");
    vbm.def_readwrite("safedist", &VBMicrolensing::corrquad,
        "Distance from planetary caustic.");


    vbm.def("LoadESPLTable", &VBMicrolensing::LoadESPLTable,
        """Loads a pre calculated binary table for extended source calculation.""");
    vbm.def("SetESPLtablefile", [](char* s) {
        VBMicrolensing::SetESPLtablefile(s);
        },
        """Sets the path to a pre calculated binary table for extended source calculation.""");
    // Maginfication calculations
    vbm.def("PSPLMag", &VBMicrolensing::PSPLMag,
        py::return_value_policy::reference,
        R"mydelimiter(
            Point Source Point Lens magnification calculation.

            Magnification of a point source by a single lens.

            Parameters
            ----------
            u : float 
                Distance of source from the center of the lens.

            Returns
            -------
            float
                Magnification.
            )mydelimiter");
    vbm.def("ESPLMag", &VBMicrolensing::ESPLMag,
        py::return_value_policy::reference,
        R"mydelimiter(
            Extended Source Point Lens magnification calculation.

            Magnification of a uniform brightness-source by a single lens.
            This uses the pre-calculated table which has to be loaded before
            calling this function.

            Parameters
            ----------
            u : float 
                Distance of source from the center of the lens.
            rho : float 
                Source radius in units of the Einstein radius of the lens.

            Returns
            -------
            float
                Magnification.
            )mydelimiter");
    vbm.def("ESPLMag2", &VBMicrolensing::ESPLMag2,
        R"mydelimiter(
            Extended Source Point Lens magnification calculation v2.0.

            ESPLMag2 works the same way as BinaryMag2. It checks whether we are
            far enough to use the point-source approximation and if necessary,
            it goes for the full computation by calling ESPLMagDark
            
            Parameters
            ----------
            u : float 
                Distance of source from the center of the lens.
            rho : float 
                Source radius in units of the Einstein radius of the lens.

            Returns
            -------
            float
                Magnification.
            )mydelimiter");
    vbm.def("ESPLMagDark", &VBMicrolensing::ESPLMagDark,
        py::return_value_policy::reference,
        R"mydelimiter(
            Extended Source Point Lens magnification calculation v2.0. 
            including limb darkening.

            Parameters
            ----------
            u : float 
                Distance of source from the center of the lens.
            rho : float 
                Source radius in units of the Einstein radius of the lens.
            a1 : float 
                Linear limb darkening coefficient. 

            Returns
            -------
            float
                Magnification.
            )mydelimiter");
    vbm.def("BinaryMag0",
        (double (VBMicrolensing::*)(double, double, double, double))
        & VBMicrolensing::BinaryMag0,
        py::return_value_policy::reference,
        R"mydelimiter(
            Magnification of a point-source by a binary lens.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 s.t. m1<m2 
            y1 : float 
                x-position of the source in the source plane.
            y2 : float 
                y-position of the source in the source plane.

            Returns
            -------
            float
                Magnification.
            )mydelimiter");

    vbm.def("BinaryMag",
        (double (VBMicrolensing::*)(double, double, double, double, double, double))
        & VBMicrolensing::BinaryMag,
        py::return_value_policy::reference,
        R"mydelimiter(
            Magnification of a uniform brightness finite source 
            by a binary lens.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 s.t. m1<m2 
            y1 : float 
                x-position of the source in the source plane.
            y2 : float 
                y-position of the source in the source plane.
            rho : float 
                Source angular radius in units of the Einstein radius 
                corresponding to the total mass.
            accuracy : float 
                Absolute accuracy goal for the magnification calculation. 

            Returns
            -------
            float
                Magnification.
            )mydelimiter");
    vbm.def("BinaryMagDark",
        &VBMicrolensing::BinaryMagDark,
        py::return_value_policy::reference,
        R"mydelimiter(
            Magnification of a limb-darkened finite source 
            by a binary lens.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 s.t. m1<m2 
            y1 : float 
                x-position of the source in the source plane.
            y2 : float 
                y-position of the source in the source plane.
            rho : float 
                Source angular radius in units of the Einstein radius 
                corresponding to the total mass.
            a1 : float 
                Source angular radius in units of the Einstein radius 
                corresponding to the total mass.

            accuracy : float 
                Absolute accuracy goal for the magnification calculation. 

            Returns
            -------
            float
                Magnification.
            )mydelimiter");
    vbm.def("BinaryMagMultiDark",
        (void (VBMicrolensing::*)(double, double, double, double, double, std::vector<double> a1_list, int, std::vector<double> mag_list, double))
        & VBMicrolensing::BinaryMagMultiDark,
        py::return_value_policy::reference,
        R"mydelimiter(
            Magnification of a limb-darkened source by a binary lens in \
            different filters with different limb darkening coefficients.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 s.t. m1<m2 
            y1 : float 
                x-position of the source in the source plane.
            y2 : float 
                y-position of the source in the source plane.
            rho : float 
                Source angular radius in units of the Einstein radius 
                corresponding to the total mass.
            a1_list : ndarray 
                Array of linear limb darkening coefficients.    
            n_filters : int 
                Number of filters. 
            mag_list : ndarray 
                Array of magnifications to be calculated by the function. 
            accuracy : float 
                Absolute accuracy goal for the magnification calculation. 

            Returns
            -------
            void
            )mydelimiter");

    vbm.def("BinaryMag2", &VBMicrolensing::BinaryMag2,
        py::return_value_policy::reference,
        R"mydelimiter(
            Magnification of a uniform brightness finite source 
            by a binary lens. New in v2.0, implements test described
            in VBMicrolensing 2.0 paper.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 s.t. m1<m2 
            y1 : float 
                x-position of the source in the source plane.
            y2 : float 
                y-position of the source in the source plane.
            rho : float 
                Source angular radius in units of the Einstein radius 
                corresponding to the total mass.
            accuracy : float 
                Absolute accuracy goal for the magnification calculation. 

            Returns
            -------
            float
                Magnification.
            )mydelimiter");

    vbm.def("ImageContours",
        [](VBMicrolensing& self, double s, double q, double y1, double y2, double rho)
        {
            _sols_for_skiplist_curve* cimages;
            std::vector<std::vector<std::vector<double> > > images;
            self.BinaryMag(s,q,y1,y2,rho,self.Tol,&cimages);
 
            for (_skiplist_curve* scurve = cimages->first; scurve; scurve = scurve->next) {
                std::vector<std::vector<double> > newimage(2);
                for (_point* scan = scurve->first; scan; scan = scan->next) {
                    newimage[0].push_back(scan->x1);
                    newimage[1].push_back(scan->x2);
                }
                images.push_back(newimage);
            }

            return images;
        },
        R"mydelimiter(
            Image contours by a binary lens.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 s.t. m1<m2 
            y1 : float 
                x-position of the source in the source plane.
            y2 : float 
                y-position of the source in the source plane.
            rho : float 
                Source angular radius in units of the Einstein radius 
                corresponding to the total mass.

            Returns
            -------
            listofimages: list[list[list[float], list[float]]]
                [images[list_of_x1,list_of_x2]]
                .
            )mydelimiter");

    vbm.def("MultiImageContours",
        [](VBMicrolensing& self, double y1, double y2, double rho)
        {
            _sols_for_skiplist_curve* cimages;
            std::vector<std::vector<std::vector<double> > > images;
            self.MultiMag(y1, y2, rho, self.Tol, &cimages);

            for (_skiplist_curve* scurve = cimages->first; scurve; scurve = scurve->next) {
                std::vector<std::vector<double> > newimage(2);
                for (_point* scan = scurve->first; scan; scan = scan->next) {
                    newimage[0].push_back(scan->x1);
                    newimage[1].push_back(scan->x2);
                }
                images.push_back(newimage);
            }

            return images;
        },
        R"mydelimiter(
            Image contours by a multiple lens. Remember to set the configuration by SetLensGeometry

            Parameters
            ----------
            y1 : float 
                x-position of the source in the source plane.
            y2 : float 
                y-position of the source in the source plane.
            rho : float 
                Source angular radius in units of the Einstein radius 
                corresponding to the total mass.

            Returns
            -------
            listofimages: list[list[list[float], list[float]]]
                [images[list_of_x1,list_of_x2]]
                .
            )mydelimiter");

    vbm.def("SetObjectCoordinates", (void (VBMicrolensing::*)(char*)) & VBMicrolensing::SetObjectCoordinates,
        R"mydelimiter(
            Sets the astronomical coordinates of the microlensing target.            
            
            Parameters
            ----------
            CoordinateString : string 
                Format \"hr:mn:sc +deg:pr:sc\".
            )mydelimiter");

    vbm.def("SetObjectCoordinates", (void (VBMicrolensing::*)(char*, char*)) & VBMicrolensing::SetObjectCoordinates,
        R"mydelimiter(
            Sets the astronomical coordinates of the microlensing target and 
            specifies the path where to look for the position tables 
            of the satellites (if any).            
            
            Parameters
            ----------
            coordinatefile : string 
                filename with astronomical coordinates.
            sattabledir : string 
                Name of the directory containing the position tables of the satellites. 
            )mydelimiter");

    // Light curve calculations
    vbm.def("PSPLLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.PSPLLightCurve(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            PSPL light curve for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_u0, log_tE, t0].
            times : list[float] 
                Array of times at which the magnification is calculated.

            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("PSPLLightCurveParallax",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.PSPLLightCurveParallax(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            PSPL light curve for a full array of observations including parallax.

            Parameters
            ----------
            params : list[float]
                List of parameters [u0, log_tE, t0, pai1, pai2] where pai1 and 
                pai2 are the components of parallax parallel and orthogonal 
                to the Earth's acceleration.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("ESPLLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.ESPLLightCurve(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            ESPL light curve for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_u0, log_tE, t0, log_rho] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");
    vbm.def("ESPLLightCurveParallax",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.ESPLLightCurveParallax(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            ESPL light curve for a full array of observations including parallax.

            Parameters
            ----------
            params : list[float]
                List of parameters [u0, log_tE, t0, log_rho, pai1, pai2] where pai1 and pai2 
                are the two components of the parallax.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");
    vbm.def("BinaryLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.BinaryLightCurve(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Static binary lens light curve for a given set of parameters.
            Uses the BinaryMag2 function.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha, log_rho, log_tE, t0]
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("BinaryLightCurveW",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.BinaryLightCurveW(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Static binary lens light curve for a given set of parameters 
            using the center of the caustic of the lens on the right as 
            a reference point for the trajectory.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0_c, alpha, log_rho, log_tE, t0_c]
                where u0_c and t0_c are defined with respect to the center of the caustic.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");
    vbm.def("BinaryLightCurveParallax",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.BinaryLightCurveParallax(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Static binary lens light curve for a given set of parameters including parallax.
            Uses the BinaryMag2 function.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha, log_rho, log_tE, t0, pai1, pai2]
                where pai1 and pai2 are the components of parallax parallel and orthogonal to the
                Earth's acceleration.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("BinaryLightCurveOrbital",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            std::vector<double> separations(times.size());
            self.BinaryLightCurveOrbital(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), separations.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s, separations };
            return results;
        },
        R"mydelimiter(
            Static binary lens light curve for a given set of parameters including parallax.
            Uses the BinaryMag2 function.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha_0, log_rho, log_tE, 
                t0, pai1, pai2, w1, w2, w3] where pai1 and pai2 are the 
                components of parallax parallel and orthogonal to the Earth's 
                acceleration. w1, w2 and w3 are the orbital parameters (assuming
                circular motion), defined as w1=(ds/dt)/s, w2=dalpha/dt, w3=(dsz/dt)/s.
            times : list[float] 
                Array of times at which the magnification is calculated.

            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array, separation-between-lenses array]
            )mydelimiter");

    vbm.def("BinaryLightCurveKepler",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            std::vector<double> separations(times.size());
            self.BinaryLightCurveKepler(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), separations.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s, separations };
            return results;
        },
        R"mydelimiter(
             binary lens light curve for a given set of parameters including keplerian orbital motion.
            Uses the BinaryMag2 function.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha_0, log_rho, log_tE, 
                t0, pai1, pai2, w1, w2, w3] where pai1 and pai2 are the 
                components of parallax parallel and orthogonal to the Earth's 
                acceleration. w1, w2 and w3 are the orbital parameters (assuming
                circular motion), defined as w1=(ds/dt)/s, w2=dalpha/dt, w3=(dsz/dt)/s.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array, separation-between-lenses array]
            )mydelimiter");

    vbm.def("BinSourceLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.BinSourceLightCurve(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Light curve for a single lens and a binary source. Sources are 
            treated as point-like.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_tE, log_fluxratio, u0_1, u0_2, t0_1, t0_2]
            times : list[float] 
                Array of times at which the magnification is calculated.
 
           Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("BinSourceLightCurveParallax",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.BinSourceLightCurveParallax(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Light curve for a single lens and a binary source including parallax.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_tE, log_fluxratio, u0_1, u0_2, t0_1, t0_2, pai1, pai2}
                where pai1 and pai2 are the components of parallax parallel and orthogonal to the
                Earth's acceleration.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("BinSourceSingleLensXallarap",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s1(times.size());
            std::vector<double> y2s1(times.size());
            std::vector<double> y1s2(times.size());
            std::vector<double> y2s2(times.size());
            self.BinSourceSingleLensXallarap(params.data(), times.data(), mags.data(),
                y1s1.data(), y2s1.data(), y1s2.data(), y2s2.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s1,y2s1,y1s2,y2s2 };
            return results;
        },
        R"mydelimiter(
            Binary source Single Lens Xallarap light curve.

            Parameters
            ----------
            params : list[float] 	
                List of parameters [log_tE, log_qs, u0, t0, xi1, xi2, 
                rho, omega, inc, phi, w2, w3] where xi1 and xi2 are the 
                components of xallarap parallel and orthogonal to the  
                seperation between the sources.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array, source 1 position y1 array, source 1 position y2 array, source 2 position y1 array,         source 2 position y2 array]
             )mydelimiter");

    vbm.def("BinSourceExtLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.BinSourceExtLightCurve(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Light curve for a single lens and a binary source. Sources are 
            treated as point-like.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_tE, log_fluxratio, u0_1, u0_2, t0_1, t0_2, rho]
            times : list[float] 
                Array of times at which the magnification is calculated.
 
           Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("BinSourceExtLightCurveXallarap",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s1(times.size());
            std::vector<double> y2s1(times.size());
            std::vector<double> y1s2(times.size());
            std::vector<double> y2s2(times.size());
            self.BinSourceExtLightCurveXallarap(params.data(), times.data(), mags.data(),
                y1s1.data(), y2s1.data(), y1s2.data(), y2s2.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s1,y2s1,y1s2,y2s2 };
            return results;
        },
        R"mydelimiter(
            Binary source light curve including xallarap for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_tE, log_FR, u01, u02, t01, t02, log_rho1, 
                                    paiN, paiE,     # components of the parallax vector
                                    w1, w2, w3,      # relative angular orbital velocity components (Einstein angle/day)
                                    ] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    source1 position y1 array, source1 position y2 array, 
                    source2 position y1 array, source2 position y2 array]
            )mydelimiter");

    vbm.def("BinSourceBinLensXallarap",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.BinSourceBinLensXallarap(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Binary source Single Lens Xallarap light curve.

            Parameters
            ----------
            params : list[float] 	
                List of parameters [log_s, log_q, u0, alpha, log_rho, log_tE, t0, xi1, xi2, 
                omega, inc, phi, log_qs] where xi1 and xi2 are the 
                components of xallarap parallel and orthogonal to the  
                separation between the sources.
            times : list[float] 
                Array of times at which the magnification is calculated.
  
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, y1 array, y2 array]
            )mydelimiter");


    vbm.def("BinSourceLightCurveXallarap",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            std::vector<double> separations(times.size());
            self.BinSourceLightCurveXallarap(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), separations.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s, separations };
            return results;
        },
        R"mydelimiter(
            Binary source light curve.

            Note that the mass ratio q between the two sources is required 
            to calculate orbital motion. Orbital motion is assumed without 
            eccentricity (see before). The parameters u0_1, u0_2, t0_1, 
            t0_2 specify the configuration at time t0 calculated as the 
            closest approach of the center of mass.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_tE, log_fluxratio, u0_1, u0_2, t0_1, t0_2, 
                pai1, pai2, q, w1, w2, w3] where pai1 and pai2 are the 
                components of parallax parallel and orthogonal to the Earth's 
                acceleration. w1, w2 and w3 are the orbital parameters (assuming
                circular motion), defined as w1=(ds/dt)/s, w2=dalpha/dt, w3=(dsz/dt)/s.
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float]] 
                [Magnification array, y1 array, y2 array, separation-between-lenses array]
            )mydelimiter");

    vbm.def("TripleLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.TripleLightCurve(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Static binary lens light curve for a given set of parameters.
            Uses the BinaryMag2 function.

            Parameters
            ----------
            params : list[float]
                List of parameters [log(s12), log(q2), u0, alpha, log(rho), log(tE), t0, log(s13), log(q3), psi]
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("TripleLightCurveParallax",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.TripleLightCurveParallax(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Static binary lens light curve for a given set of parameters.
            Uses the BinaryMag2 function.

            Parameters
            ----------
            params : list[float]
                List of parameters [log(s12), log(q2), u0, alpha, log(rho), log(tE), t0, log(s13), log(q3), psi,px1,px2]
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("LightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.LightCurve(params.data(), times.data(), mags.data(),
                y1s.data(), y2s.data(), times.size(), (params.size() - 4) / 3 + 1);
            std::vector< std::vector<double> > results{ mags,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Static binary lens light curve for a given set of parameters.
            Uses the BinaryMag2 function.

            Parameters
            ----------
            params : list[float]
                List of parameters [t0, log_tE, log_rho, s1_im, s2_real,....,s2_im,...., q2,...,qn]
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float]] 
                [Magnification array, source position y1 array, source position y2 array]
            )mydelimiter");

    // Astrometric functions

    vbm.def("CombineCentroids",
        [](VBMicrolensing& self, std::vector< std::vector<double> >  results, double blending)
        {
            std::vector<double> c1tot(results[0].size());
            std::vector<double> c2tot(results[0].size());
            self.CombineCentroids(results[0].data(), results[1].data(), results[2].data(), results[3].data(), results[4].data(), c1tot.data(), c2tot.data(), blending, results[0].size());
            std::vector< std::vector<double> > centroids{ c1tot, c2tot };
            return centroids;
        },
        R"mydelimiter(
            Combine source and lens centroids from astrometric functions taking into account blending and magnification.

            Parameters
            ----------
            results : list[list[float],list[float],list[float],list[float],list[float],list[float],list[float]]
                Results returned from an astrometric functions.
            blending : float 
                flux ratio lens/source.
 
            Returns
            -------
            centroids: list[list[float],list[float]] 
                [combined centroid of images Dec array, combined centroid of images R.A. array]
            )mydelimiter");

    vbm.def("PSPLAstroLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> c1s(times.size());
            std::vector<double> c2s(times.size());
            std::vector<double> c1l(times.size());
            std::vector<double> c2l(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.astrometry = true;
            self.PSPLAstroLightCurve(params.data(), times.data(), mags.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags, c1s, c2s, c1l, c2l,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            PSPL light curve and astrometry for a full array of observations.

            Parameters
            ----------
            params : list[float]
                 List of parameters [u0, log_tE, t0, 
                                    paiN, paiE,     #components of the parallax vector
                                    muS_N, muS_E,   # proper motion components of the source (mas/yr)
                                    pai_S,          # parallax of the source (mas)
                                    thetaE          # Einstein angle (mas) 
             times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    centroid of images N array, centroid of images E array, 
                    centroid of lens N array, centroid of lens E array,
                    source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("ESPLAstroLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> c1s(times.size());
            std::vector<double> c2s(times.size());
            std::vector<double> c1l(times.size());
            std::vector<double> c2l(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.astrometry = true;
            self.parallaxsystem = 1;
            self.ESPLAstroLightCurve(params.data(), times.data(), mags.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags, c1s, c2s, c1l, c2l,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            ESPL light curve and astrometry for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [u0, log_tE, t0, log_rho, 
                                    paiN, paiE,     #components of the parallax vector
                                    muS_N, muS_E,   # proper motion components of the source (mas/yr)
                                    pai_S,          # parallax of the source (mas)
                                    thetaE          # Einstein angle (mas) 
                                    ] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    centroid of images N array, centroid of images E array, 
                    centroid of lens N array, centroid of lens E array,
                    source position y1 array, source position y2 array]
            )mydelimiter");


    vbm.def("BinaryAstroLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> c1s(times.size());
            std::vector<double> c2s(times.size());
            std::vector<double> c1l(times.size());
            std::vector<double> c2l(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.astrometry = true;
            self.parallaxsystem = 1;
            self.BinaryAstroLightCurve(params.data(), times.data(), mags.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags, c1s, c2s, c1l, c2l,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Binary light curve and astrometry for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha, log_rho, log_tE, t0, 
                                    paiN, paiE,     #components of the parallax vector
                                    muS_N, muS_E,   # proper motion components of the source (mas/yr)
                                    pai_S,          # parallax of the source (mas)
                                    thetaE          # Einstein angle (mas) 
                                    ] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    centroid of images N array, centroid of images E array, 
                    centroid of lens N array, centroid of lens E array,
                    source position y1 array, source position y2 array]
            )mydelimiter");

    vbm.def("BinaryAstroLightCurveOrbital",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> c1s(times.size());
            std::vector<double> c2s(times.size());
            std::vector<double> c1l(times.size());
            std::vector<double> c2l(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            std::vector<double> seps(times.size());
            self.astrometry = true;
            self.parallaxsystem = 1;
            self.BinaryAstroLightCurveOrbital(params.data(), times.data(), mags.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                y1s.data(), y2s.data(), seps.data(), times.size());
            std::vector< std::vector<double> > results{ mags, c1s, c2s, c1l, c2l,y1s,y2s, seps };
            return results;
        },
        R"mydelimiter(
            Binary light curve and astrometry including circular orbital motion for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha, log_rho, log_tE, t0, 
                                    paiN, paiE,     # components of the parallax vector
                                    w1, w2, w3,      # relative angular orbital velocity components (Einstein angle/day)
                                    muS_N, muS_E,   # proper motion components of the source (mas/yr)
                                    pai_S,          # parallax of the source (mas)
                                    thetaE          # Einstein angle (mas) 
                                    ] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    centroid of images N array, centroid of images E array, 
                    centroid of lens N array, centroid of lens E array,
                    source position y1 array, source position y2 array, separations between the lenses array]
            )mydelimiter");

    vbm.def("BinaryAstroLightCurveKepler",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> c1s(times.size());
            std::vector<double> c2s(times.size());
            std::vector<double> c1l(times.size());
            std::vector<double> c2l(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            std::vector<double> seps(times.size());
            self.astrometry = true;
            self.parallaxsystem = 1;
            self.BinaryAstroLightCurveKepler(params.data(), times.data(), mags.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                y1s.data(), y2s.data(), seps.data(), times.size());
            std::vector< std::vector<double> > results{ mags, c1s, c2s, c1l, c2l,y1s,y2s, seps };
            return results;
        },
        R"mydelimiter(
            Binary light curve and astrometry including eccentric orbital motion for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha, log_rho, log_tE, t0, 
                                    paiN, paiE,     # components of the parallax vector
                                    w1, w2, w3,      # relative angular orbital velocity components (Einstein angle/day)
                                    sz_s,          # Ratio of separation along the line of sight and the transverse separation at t0
                                    a_stot,         # Semimajor axis over the 3D separation at time t0
                                    muS_N, muS_E,   # proper motion components of the source (mas/yr)
                                    pai_S,          # parallax of the source (mas)
                                    thetaE          # Einstein angle (mas) 
                                    ] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    centroid of images N array, centroid of images E array, 
                    centroid of lens N array, centroid of lens E array,
                    source position y1 array, source position y2 array, separations between the lenses array]
            )mydelimiter");

    vbm.def("BinSourceAstroLightCurveXallarap",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> c1s(times.size());
            std::vector<double> c2s(times.size());
            std::vector<double> c1l(times.size());
            std::vector<double> c2l(times.size());
            std::vector<double> y1s1(times.size());
            std::vector<double> y2s1(times.size());
            std::vector<double> y1s2(times.size());
            std::vector<double> y2s2(times.size());
            self.astrometry = true;
            self.parallaxsystem = 1;
            self.BinSourceAstroLightCurveXallarap(params.data(), times.data(), mags.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                y1s1.data(), y2s1.data(), y1s2.data(), y2s2.data(), times.size());
            std::vector< std::vector<double> > results{ mags, c1s, c2s, c1l, c2l,y1s1,y2s1,y1s2,y2s2 };
            return results;
        },
        R"mydelimiter(
            Binary source light curve and astrometry including xallarap for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_tE, log_FR, u01, u02, t01, t02, log_rho1, 
                                    paiN, paiE,     # components of the parallax vector
                                    w1, w2, w3,      # relative angular orbital velocity components (Einstein angle/day)
                                    muS_N, muS_E,   # proper motion components of the source (mas/yr)
                                    pai_S,          # parallax of the source (mas)
                                    thetaE          # Einstein angle (mas) 
                                    ] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    centroid of images N array, centroid of images E array, 
                    centroid of lens N array, centroid of lens E array,
                    source1 position y1 array, source1 position y2 array, 
                    source2 position y1 array, source2 position y2 array]
            )mydelimiter");


    vbm.def("TripleAstroLightCurve",
        [](VBMicrolensing& self, std::vector<double> params, std::vector<double> times)
        {
            std::vector<double> mags(times.size());
            std::vector<double> c1s(times.size());
            std::vector<double> c2s(times.size());
            std::vector<double> c1l(times.size());
            std::vector<double> c2l(times.size());
            std::vector<double> y1s(times.size());
            std::vector<double> y2s(times.size());
            self.astrometry = true;
            self.parallaxsystem = 1;
            self.TripleAstroLightCurve(params.data(), times.data(), mags.data(), c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                y1s.data(), y2s.data(), times.size());
            std::vector< std::vector<double> > results{ mags, c1s, c2s, c1l, c2l,y1s,y2s };
            return results;
        },
        R"mydelimiter(
            Triple light curve and astrometry for a full array of observations.

            Parameters
            ----------
            params : list[float]
                List of parameters [log_s, log_q, u0, alpha, log_rho, log_tE, t0, 
                                    log(s13), log(q3), psi
                                    paiN, paiE,     #components of the parallax vector
                                    muS_N, muS_E,   # proper motion components of the source (mas/yr)
                                    pai_S,          # parallax of the source (mas)
                                    thetaE          # Einstein angle (mas) 
                                    ] 
            times : list[float] 
                Array of times at which the magnification is calculated.
 
            Returns
            -------
            results: list[list[float],list[float],list[float],list[float],list[float],list[float],list[float]] 
                [Magnification array,
                    centroid of images N array, centroid of images E array, 
                    centroid of lens N array, centroid of lens E array,
                    source position y1 array, source position y2 array]
            )mydelimiter");
    // Other functions


    vbm.def("Multicaustics",
        [](VBMicrolensing& self)
        {
            _sols* critcau;

            critcau = self.PlotCrit();
            int ncaus = critcau->length / 2;
            std::list <std::vector<std::vector<double>>> caustics{};
            _curve* c = critcau->first;
            for (int i = 0; i < ncaus; i++) c = c->next;
            for (int i = 0; i < ncaus; i++) {
                std::vector<double> y(c->length);
                std::vector<std::vector<double>> cau(2, y);
                _point* p = c->first;
                for (int j = 0; j < c->length; j++) {
                    cau[0][j] = p->x1;
                    cau[1][j] = p->x2;
                    p = p->next;
                }
                caustics.push_back(cau);
                c = c->next;
            }
            delete critcau;
            return caustics;
        },
        R"mydelimiter(
            Caustics for given separation and mass ratio.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 such that m1<m2 

            Returns
            -------
            solutions : _sols
                List of caustics.
            )mydelimiter");

    vbm.def("Multicriticalcurves",
        [](VBMicrolensing& self)
        {
            _sols* critcau;

            critcau = self.PlotCrit();
            int ncrits = critcau->length / 2;
            std::list <std::vector<std::vector<double>>> criticalcurves{};
            _curve* c = critcau->first;
            for (int i = 0; i < ncrits; i++) {
                std::vector<double> y(c->length);
                std::vector<std::vector<double>> crit(2, y);
                _point* p = c->first;
                for (int j = 0; j < c->length; j++) {
                    crit[0][j] = p->x1;
                    crit[1][j] = p->x2;
                    p = p->next;
                }
                criticalcurves.push_back(crit);
                c = c->next;
            }
            delete critcau;
            return criticalcurves;
        },
        R"mydelimiter(
            Critical curves for given separation and mass ratio.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 such that m1<m2 

            Returns
            -------
            solutions : _sols
                List of critical curves.
            )mydelimiter");
    vbm.def("Caustics",
        [](VBMicrolensing& self, double s, double q)
        {
            _sols* critcau;

            critcau = self.PlotCrit(s, q);
            int ncaus = critcau->length / 2;
            std::list <std::vector<std::vector<double>>> caustics{};
            _curve* c = critcau->first;
            for (int i = 0; i < ncaus; i++) c = c->next;
            for (int i = 0; i < ncaus; i++) {
                std::vector<double> y(c->length);
                std::vector<std::vector<double>> cau(2, y);
                _point* p = c->first;
                for (int j = 0; j < c->length; j++) {
                    cau[0][j] = p->x1;
                    cau[1][j] = p->x2;
                    p = p->next;
                }
                caustics.push_back(cau);
                c = c->next;
            }
            delete critcau;
            return caustics;
        },
        R"mydelimiter(
            Caustics for given separation and mass ratio.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 such that m1<m2 

            Returns
            -------
            solutions : _sols
                List of caustics.
            )mydelimiter");

    vbm.def("Criticalcurves",
        [](VBMicrolensing& self, double s, double q)
        {
            _sols* critcau;

            critcau = self.PlotCrit(s, q);
            int ncrits = critcau->length / 2;
            std::list <std::vector<std::vector<double>>> criticalcurves{};
            _curve* c = critcau->first;
            for (int i = 0; i < ncrits; i++) {
                std::vector<double> y(c->length);
                std::vector<std::vector<double>> crit(2, y);
                _point* p = c->first;
                for (int j = 0; j < c->length; j++) {
                    crit[0][j] = p->x1;
                    crit[1][j] = p->x2;
                    p = p->next;
                }
                criticalcurves.push_back(crit);
                c = c->next;
            }
            delete critcau;
            return criticalcurves;
        },
        R"mydelimiter(
            Critical curves for given separation and mass ratio.

            Parameters
            ----------
            s : float 
                The projected separation of the binary lens in units of the 
                Einstein radius corresponding to the total mass.
            q : float 
                Binary lens mass fraction q = m1/m2 such that m1<m2 

            Returns
            -------
            solutions : _sols
                List of critical curves.
            )mydelimiter");

    // Limb darkening
    /*vbm.def("SetLDprofile", (void (VBMicrolensing::*)(double(*)(double), int))
            &VBMicrolensing::SetLDprofile,
            "Set LD profile using a user-defined function");*/

    vbm.def("SetLDprofile", (void (VBMicrolensing::*)(VBMicrolensing::LDprofiles))
        & VBMicrolensing::SetLDprofile,
        "Set LD profile using a predefined profile");

    vbm.def("SetLensGeometry",
        [](VBMicrolensing& self, std::vector<double> pr) {
            self.SetLensGeometry(pr.size() / 3, pr.data());
        },
        py::return_value_policy::reference,
        R"mydelimiter(
            Set the geometry of the system

            Parameters
            ----------
            pr : list[float]
                List of parameters [z1_re, z1_im, q1, z2_re, z2_im, q2,....,zn_re,zn_im,qn]
                
            )mydelimiter");


    vbm.def("MultiMag0",
        [](VBMicrolensing& self, double y1, double y2) -> double {
            return self.MultiMag0(y1, y2);
        },
        py::return_value_policy::reference,
        R"mydelimiter(
            Magnification of a point-source by a multiple lens.

            Parameters
            ----------
            source : tuple of floats
                A tuple (y1, y2) where y1 is the real-position and y2 is the imaginary-position 
                of the source in the source plane.

            Returns
            -------
            float
                Magnification.
            )mydelimiter");

    vbm.def("MultiMag",
        [](VBMicrolensing& self, double y1, double y2, double rho) -> double {
            return self.MultiMag(y1, y2, rho);
        },
        py::return_value_policy::reference,
        R"pbdoc(
        Compute the magnification of a uniform brightness finite source 
        by a multiple lens.

        Parameters
        ----------
        y1 : float 
                x-position of the source in the source plane.
        y2 : float 
            y-position of the source in the source plane.
        rho : float 
            The source angular radius in units of the Einstein radius 
            corresponding to the total mass.
        

        Returns
        -------
        float
            The magnification.
        )pbdoc");


    vbm.def("MultiMagDark",
        [](VBMicrolensing& self, double y1, double y2, double rho, double Tol) -> double {
            return self.MultiMagDark(y1, y2, rho, Tol);
        },
        py::return_value_policy::reference,
        R"pbdoc(
        Magnification of a limb-darkened finite source 
        by a multiple lens.

        Parameters
        ----------
        y1 : float 
                x-position of the source in the source plane.
        y2 : float 
            y-position of the source in the source plane.
        rho : float 
            The source angular radius in units of the Einstein radius 
            corresponding to the total mass.
        accuracy : float 
                Absolute accuracy goal for the magnification calculation.

        Returns
        -------
        float
            The magnification.
        )pbdoc");

    vbm.def("MultiMag2",
        [](VBMicrolensing& self, double y1, double y2, double rho) -> double {
            return self.MultiMag2(y1, y2, rho);
        },
        py::return_value_policy::reference,
        R"pbdoc(
        Compute the magnification of a finite source 
        by a multiple lens. In v2.0, implements test described
            in VBMicrolensing 2.0 paper.

        Parameters
        ----------
        y1 : float 
                x-position of the source in the source plane.
        y2 : float 
            y-position of the source in the source plane.
        rho : float 
            The source angular radius in units of the Einstein radius 
            corresponding to the total mass.
        

        Returns
        -------
        float
            The magnification.
        )pbdoc");

    vbm.def("SetMethod",
        &VBMicrolensing::SetMethod,
        "User choice of Method");

    //  Method: Singlepoly, Multipoly, Nopoly
    py::enum_<VBMicrolensing::Method>(vbm, "Method")
        .value("Singlepoly", VBMicrolensing::Method::Singlepoly)
        .value("Multipoly", VBMicrolensing::Method::Multipoly)
        .value("Nopoly", VBMicrolensing::Method::Nopoly)
        .export_values();

    //LDlinear, LDquadratic, LDsquareroot, LDlog, LDuser
    py::enum_<VBMicrolensing::LDprofiles>(vbm, "LDprofiles")
        .value("LDlinear", VBMicrolensing::LDprofiles::LDlinear)
        .value("LDquadratic", VBMicrolensing::LDprofiles::LDquadratic)
        .value("LDsquareroot", VBMicrolensing::LDprofiles::LDsquareroot)
        .value("LDlog", VBMicrolensing::LDprofiles::LDlog)
        .value("LDuser", VBMicrolensing::LDprofiles::LDuser)
        .export_values();

    py::class_<_theta>(m, "_theta")
        .def(py::init<double>()); //constructor 

    py::class_<_point>(m, "_point")
        .def(py::init<double, double, _theta*>())
        .def_readwrite("next", &_point::next)
        .def_readwrite("prev", &_point::prev)
        .def_readonly("x1", &_point::x1)
        .def_readonly("x2", &_point::x2);

    py::class_<_curve>(m, "_curve")
        .def(py::init<_point*>()) //constructor 1
        .def(py::init()) //constructor 2
        .def_readwrite("first", &_curve::first)
        .def_readwrite("last", &_curve::last)
        .def_readwrite("next", &_curve::next)
        .def_readwrite("prev", &_curve::prev);

    py::class_<_skiplist_curve>(m, "_skiplist_curve")
        .def(py::init<_point*,int>()) //constructor 1
        .def(py::init()) //constructor 2
        .def_readwrite("first", &_skiplist_curve::first)
        .def_readwrite("last", &_skiplist_curve::last)
        .def_readwrite("next", &_skiplist_curve::next)
        .def_readwrite("prev", &_skiplist_curve::prev);

    py::class_<_sols>(m, "_sols")
        .def(py::init()) //constructor
        .def_readwrite("first", &_sols::first)
        .def_readwrite("last", &_sols::last);

    py::class_<_sols_for_skiplist_curve>(m, "_sols_for_skiplist_curve")
        .def(py::init()) //constructor
        .def_readwrite("first", &_sols_for_skiplist_curve::first)
        .def_readwrite("last", &_sols_for_skiplist_curve::last);

    py::class_<complex>(m, "complex")
        .def(py::init<double, double>())
        .def(py::init<double>())
        .def(py::init<>())
        .def_readwrite("re", &complex::re)
        .def_readwrite("im", &complex::im);
}
