# GalSim
### Generates a population of stars in the Milky Way.

******************************

####Overview
GalSim is a Mily Way stellar population generator. Using the electron density distribution of [Cordes & Lazio (2002)](#C1), the Schmidt-Kennicut law ([Kennicut 1998](#K98)), the initial mass function of [Kroupa (2001)](#K1), and the accretion model by [McKee & Tan (2003)](#M3) this code produces a distribution of stars each with an age, current stellar mass and coordinates within the Galaxy.

####Example Usage
After building, the directory tree (with `bin` as root directory) of the application should look like this:
```
.
├── galsim
├── config
|   ├── required
|   |   ├── hosokawa_tracks_interp
|   |   |   ├── int_md.dat
|   |   |   └── mass_ms_fine.dat
|   |   ├── cepheids.csv
|   |   ├── cepheids_pos.csv
|   |   ├── clust_ak.txt
|   |   ├── lyman_101108_corr-highmass.dat
|   |   └── msx_data.csv
|   └── usr
|       └── galsim-config.lua
```

To run, execute `galsim`:
```bash
./galsim
```
By default GalSim reads in the configuration file, `config/usr/galsim-config.lua`.
You can specify your own configuration file:
```bash
./galsim --paramfile=/path/to/galsim-config.lua
```

#####Setup

For example, parameters could be:

```lua
-- config/usr/galsim-config.lua

Parameters = {
	densityfile =                "",
	starpopfile =                "",
	output_directory =           "tmp",

	sfr =                        1.5,
	total_time =                 1.0e6,

	radius =                     18.0,
	height =                     2.0,
	resolution =                 0.1,
	solar_position_x =           0.0,
	solar_position_y =           8.5,
	solar_position_z =           0.0,
	A_inner =                    5.0,
	A_a =                        11.0,
	A_1 =                        17.0,
	A_2 =                        3.7,
	H_1 =                        0.95,
	H_2 =                        0.14,
	n_1 =                        0.0347,
	n_2 =                        0.09,
	h_1 =                        0.95,

	wa_w1 =                      0.60,
	wa_w2 =                      0.90,
	wa_w3 =                      0.60,
	wa_w4 =                      0.48,
	wa_w5 =                      0.60,

	na_f1 =                      0.0150,
	na_f2 =                      0.0360,
	na_f3 =                      0.0390,
	na_f4 =                      0.0300,
	na_f5 =                      0.0075,

	ha_h1 =                      0.250,
	ha_h2 =                      0.200,
	ha_h3 =                      0.325,
	ha_h4 =                      0.375,
	ha_h5 =                      0.250,

	arm_joins_1 =                45.0,
	arm_joins_2 =                45.0,
	arm_joins_3 =                -135.0,
	arm_joins_4 =                -135.0,

	arm_min_radius =             1.0,

	logarm_rmin_1 =              4.0,
	logarm_rmin_2 =              4.0,
	logarm_rmin_3 =              4.0,
	logarm_rmin_4 =              4.0,

	logarm_thmin_1 =             131.761,
	logarm_thmin_2 =             37.214,
	logarm_thmin_3 =             -61.889,
	logarm_thmin_4 =             -139.832,

	logarm_asp_1 =               524.677,
	logarm_asp_2 =               482.624,
	logarm_asp_3 =               572.652,
	logarm_asp_4 =               551.500,

	logarm_rmax_1 =              13.0,
	logarm_rmax_2 =              13.0,
	logarm_rmax_3 =              13.0,
	logarm_rmax_4 =              13.0,
}
```

#####Output
GalSim outputs data files to the directory assigned to `output_directory`.

`density3D.txt` is the Galactic number density grid data. The three integers on the first row are the resolutions of the Galactic grid along the x, y, and z axes respectively. On the third row starts the column data. From left to right the data are: Galactic x, y, and z coordinates; thin, thick, and spiral arm number density components; and finally the total number density in that grid cell.

`density-xy.txt` is the Galactic number densities summed along the z-axis. The first two integers on the first row are the resolutions of the Galactic grid along the x, y, and z axes respectively. On the third row stars the column data. From left to right the data are: Galactic x and y coordinates; thin, thick, and spiral arm number density components; and finally the total summed number density in that grid cell.

`density-xz.txt` is the same as `density-xy.txt` except the data is summed along the y-axis.

`starpop.txt` is the star data before we have utilised an accretion model. From left to right the columns are: final stellar mass; age; distance from the sun; distance from the Galactic centre; column density; V-band, K-band and 21um extinctions; Galactic x, y, and z coordinates; and Galactic longitude and latitude.

`starpopfinal.txt` is the star data after we have utilised an accretion model. From left to right the columns are: current and final stellar masses; time on main-sequence; age; distance from the sun; distance from the Galactic centre; column density; V-band, K-band and 21um extinctions; bolometric luminosity; 21um flux level; Galactic x, y, and z coordinates; and Galactic longitude and latitude.  

####Compiling

GalSim uses the cmake build process. To build simply make a `build` directory and call `ccmake` from there:
```bash
mkdir build
cd mkdir
ccmake path/to/GalSim
make
```

To specify your C++ compiler (GalSim requires gcc 4.7.2+) you need to set some environment variables before calling `ccmake`:
```bash
export CC=/path/to/gcc
export CXX=/path/to/g++
```

####Advanced Usage
The parameters not included in this table should not be modified unless you know what you're doing.

#####Basic
| Parameter                     | Notes                                     |
| :---------------------------- | :---------------------------------------- |
| ```densityfile```             | 3D density data that is output by this program (so the program doesn't have to recalculate it). |
| ```starpopfile```             | Stellar population file. For recalculating stellar masses predicted by another accretion model. |
| ```output_directory```        | The output directory. |
| ```sfr```                     | Star formation rate in units of M_solar per year. |
| ```total_time```              | Maximum age of a star in the distribution. |

####Developer info
Harrison Steggles, University of Leeds (PhD student).

####References
<a name="C1"></a>Cordes, J. M. and Lazio, T. J. W. (2002). NE2001.I. A New Model for the Galactic Distribution of Free Electrons and its Fluctuations. ArXiv Astrophysics e-prints.([link](http://arxiv.org/pdf/astro-ph/0207156v3.pdf))  
<a name="K98"></a>Kennicutt, Jr., R. C. (1998). The Global Schmidt Law in Star-forming Galaxies. ApJ, 498:541–552.([link](http://iopscience.iop.org/article/10.1086/305588/pdf))  
<a name="K1"></a>Kroupa, P. (2001). On the variation of the initial mass function. MNRAS, 322:231–246.([link](http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?2001MNRAS.322..231K&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf))  
<a name="M3"></a>McKee, C. F. and Tan, J. C. (2003). The Formation of Massive Stars from Turbulent Cores. ApJ, 585:850–871.([link](http://iopscience.iop.org/article/10.1086/346149/pdf))  


####Requirements
* gcc 4.7.2+.
* [Selene](https://github.com/jeremyong/Selene): "Simple C++11 friendly header-only bindings to Lua 5.2+".  
* [Lua5.2+](http://www.lua.org/): "A powerful, fast, lightweight, embeddable scripting language".  
