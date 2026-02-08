Sunangle
========

This crate calculates the position of the sun in the sky, and other
related things like sunrise time.


Accuracy
========

This is a new crate, and also I could mess it up at any time.

But at some point the data was checked against spreadsheet calculations
from [NOAA](https://gml.noaa.gov/grad/solcalc/calcdetails.html)
which use formulas published in `Astronomical Algorithms` by
Jean Meeus.

Some alternate newer formulas may be more accurate, and some are
available (and more will become available).

I have not yet coded enough alternate formulas to match and verify
data against popular sources such as https://www.suncalc.org


Is it a complicated formula?
============================

That is putting it mildly.

First, you need to know the exact time.  Not in your timezone, but
at GMT.  And not with any daylight savings time adjustments. Oh, but
also not in UTC even at all but rather in Tt (terresterial time)
which 32.184 seconds plus 9 plus a number of additional leap seconds
ahead of UTC.  For this, we use another crate of mine
[astrotime](https://github.com/mikedilger/astrotime).

Then there are a whole bunch of calculations that go something like
this (details elided, read the code):

* Julian Day
* Julian Century in J2000
* Solar geometric mean longitude
* Solar geometric mean anomoly
* Earth orbit eccentricitiy
* Equasion of center
* True longitude
* True anomoly
* Rad vector in AU (optional)
* Apparent longitude
* Obliquity
* Right ascension
* Declination
* Equasion of Time
  * Sunrise (optional)
  * Solar Noon in UTC (optional)
  * Sunrise in UTC (optional)
  * Sunset in UTC (optional)
  * Sunlight duration minutes (optional)
* True solar time in minutes
* Hour angle
* Zenith angle
* Elevation angle
* Approximate atmospheric refraction
* Corrected zenith angle
* Corrected elevation angle
* Azimuth angle
* Vector
