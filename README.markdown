Chronos
=======

> I have witnessed Your beauty with my eyes and my mind and I have witnessed Your complexity with numbers.

This library is intended for astronomical purposes and provides functionality realted to observations of the sky.

Some of the functions offered by this library compute:

* **Julian date** and **Julian Ephemeris date**
* Difference between **Dynamical Time** and **Universal Time**, commonly denoted as **ΔT**
* Date of **Easter**
* Conversion between **horizontal**, **equatorial** and **ecliptic** coordinates
* Times of **rising** and **setting** of celesetial bodies
* Corrections due to effects of **precession** and **nutation** of the Earth's rotational axis, **aberration** of light,
  diurnal **parallax** and atmospheric **refraction**
* Times of **solstices** and **equinoxes**
* **True** and **apparent positions** of the Sun, Moon and major planets of the Solar System
* Apparent **illuminated fractions** of the disks of the Moon and the planets and their **apparent magnitudes**

![Chronos logo][0]

Refer to the comments in each header file for more information on full functionality of the library as well as
information on the functions themselves.

The main source material during developement was an amazing book by Belgian astronomer and scientist J. Meeus called
["Astronomical algorithms"][1]. Howerver, whenever this book referenced another material, use of original articles was
made. For instance, original paper of 1980 IAU Theory of Nutation was used when writing code of correction for nutation
of the rotational axis of the Earth, as well as original papers by L.V. Morrison and F.R. Stephenson when writing code
that computes the difference between Dynamical Time and Universal Time, and some other.

Code is thoroughly documented. Header file comments focus on explanations of how to use provided functions, what
parameters are to be substituded, in what units, reference frames, etc. Each header file contains a brief comment on
top which briefly introduces terms related to corrsponding routines. Programmers with very basic knowledge of astronomy
should feel very comfortable with provided documentation.

Every function references source material that was used during its implementation.

Source files contain comments mainly related to the algorithms and methodologies themselves and assume a proper
knowledge of astronomy to some extent. This should help making improvements to the library much easier as well as
maximizing the ability of others to understand algorithms that have been put into these functions.

One of the focuses was on code simplicty and readability. In some places it may feel bloated or unnecesary, but the
idea was to preserve the steps involved in algorithms and make them explicit.

Library uses planetary theory VSOP87 version D and lunar theory ELP version ELP2000-82B to compute positions of the
planets (and the Sun) and the Moon respectively. This theories have been implemented separately as standalone libraries
which Chronos must be linked agains. Source code for the libraries is avalable at these links:

* [VSOP87D][2]
* [ELP2000-82B][3]

At the given links you'll find more information about these theories as well as links to original papers by
P. Bretagnon and G. Francou for VSOP and M. Chapront-Touzé and J. Chapront for ELP.

<br/>
**If you find this library useful and would like to see additional functionality implemented, your offers are very
welcome. Also, if you find any errors, mistakes, typos or inconsistencies, please, contact me via [e-mail][4].**

<br />
Serhii Tsyba

18.08.2010<br />
Helsinki, Finland

---

Staged future improvements and expected functionality
-----------------------------------------------------

* Move from ELP theory version ELP2000-82B to version ELP/MPP02 for calculations of lunar positions
* Move from VSOP theory version VSOP87 to VSOP2000 for calculations of planetary positions and positions of the Sun (?)
* Move from IAU 1980 Theory of Nutation to IAU 2000 Theory of Nutation
* Move nutation code into a separate library (?)
* Move code that computes mean orbital elements to the VSOP library
* Add calculations of positions of minor planets (elliptic motion)
* Add calculations of positions of comets (parabolic motion)
* Add calculations of times of solar and lunar eclipses
* Improve calculations of rising and setting times
* Add some related calculations, like solving equation of Kepler, solving equation of the centre, etc.


[0]: https://dl.dropbox.com/u/4936034/Referred/chronos.jpg      "Chronos logo"
[1]: http://www.willbell.com/math/mc1.HTM                       "J. Meeus. Astronomical algorithms"
[2]: http://github.com/sertsy/vsop87d                           "VSOP87D sources at GitHub"
[3]: http://github.com/sertsy/elp2000-82b                       "ELP2000-82B sources at GitHub"
[4]: mailto:sertsy@gmail.com                                    "Contact e-mail"