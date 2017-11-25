# Programmatic solution to question from Quora.com about minute and hour hands on a clock.

The original question is:

How many times are the minute hand and hour hands are [sic] at 30 degrees in a day?

The Python script, min_hour_30deg.py, uses the [SpiceyPy package](https://github.com/AndrewAnnex/SpiceyPy) to simulate a clock using three synthetic bodies:  a central body; two bodies orbiting the central body.  The central body represents the clock; the two bodies represent the minute and hour hands.  The trajectories of the minute and hour bodies are circular orbits around the clock body such that

* the minute body completes each orbit of the clock once per minute,
* the hour body completes each orbit of the clock once per hour,
* the minute and hour bodies are at a Right Ascension of zero degrees at 2001-JAN-01-00:00:00 TDB (TDB => Barycentric Dynamical Time).

## Manifest

### min_hour_30deg.py

* Script to verify there are 44 occurences per day of times when hour and minute hand are 30degrees apart.
* Usage:  python min_hour_30deg.py [--debug]
* This script uses the [SpiceyPy package](https://github.com/AndrewAnnex/SpiceyPy) to model the MINUTE and HOUR hands as bodies orbiting a CLOCK body; the angle between the minute and hour hands are mmodeled as the [phase angle](https://en.wikipedia.org/wiki/Phase_angle_(astronomy)) with the CLOCK as the target body, the MINUTE body as the illuminating body, and the HOUR body as the observing body; the phase angle is the angle between the [CLOCK=>MINUTE] and [CLOCK=>HOUR] vectors.

### test_output.txt

* Output of script, with debugging turned on

### min_hour_clock.bsp

* [NAIF/JPL SPICE](http://naif.jpl.nasa.gov) SP-Kernel containing synthetic trajectories of MINUTE and HOUR bodies that orbit CLOCK body
* Script will recreate this file in the current working directory if it does not exist

### README.md

* This file
