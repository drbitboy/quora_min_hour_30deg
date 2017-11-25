# Programmatic solution to question from Quora.com about minut and hour hands on a clock.

The original question is:

How many times are the minute hand and hour hands are [sic] at 30 degrees in a day?

The Python script, min_hour_30deg.py, uses the [SpiceyPy package](https://github.com/AndrewAnnex/SpiceyPy) to simulate a clock using three synthetic bodies:  a central body; two bodies orbiting the central body.  The central body represents the clock; the two bodies represent the minute and hour hands.  The trajectories of the minute and hour bodies are circular orbits around the clock body such that

* the minute body completes each orbit of the clock once per minute,
* the hour body completes each orbit of the clock once per hour,
* the minute and hour bodies are at a Right Ascension of zero degrees at 2001-JAN-01-00:00:00 TDB (TDB => Barycentric Dynamical Time).
