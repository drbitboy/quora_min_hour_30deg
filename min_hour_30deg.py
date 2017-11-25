"""
\begindata
  NAIF_BODY_CODE += ( 2000000,  2000060, 2003600 )
  NAIF_BODY_NAME += ( 'CLOCK', 'MINUTE',  'HOUR' )

  BODY2000000_POLE_RA  = ( 270 0 0 )
  BODY2000000_POLE_DEC = ( 90. 0 0 )
  BODY2000000_POLE_PM  = ( 0. 0 0 )

  ET0 = @2000-01-01-00:00:00.000

  CLOCKSPK = 'min_hour_clock.bsp'
\begintext
"""
import os
import sys
import spiceypy as sp

if "__main__" == __name__:

  ### FURNSH (SPICE furnish i.e. load) this script as a kernel, as well as
  ### any command-line arguments
  map(sp.furnsh,[arg for arg in sys.argv if arg != '--debug'])
  doDebug = '--debug' in sys.argv

  ### SPICE names and IDs of Clock, minute hand, hour hand
  sClock, sMinute, sHour = sNames = "CLOCK MINUTE HOUR".split()
  iClock, iMinute, iHour = iIds = map(sp.bods2c,sNames)

  ### SPICE works in seconds; get seconds per minute, per hour, per day
  spm = sp.convrt(1.0, 'minutes', 'seconds')
  sph = sp.convrt(1.0, 'hours', 'seconds')
  hpd = sp.convrt(1.0, 'days', 'hours')
  spd = sp.spd()

  ### Other constants:  2*PI; degrees/radian
  twopi = sp.twopi()
  dpr = sp.dpr()

  ### Get start time (TDB) and SPK filename from kernel pool
  et0 = sp.gdpool('ET0',0,1)[0]
  fnClockSpk = sp.gcpool('CLOCKSPK', 0, 1, 99)[0]

  ### Make clock SPK if it does not exist
  if not os.path.exists(fnClockSpk):

    ### Orbital period formula:
    ###   T = 2 PI (a*3 / mu))**(1/2)
    ### where
    ###    T = orbital period
    ###    a = length of semi-major axis
    ###   mu = GM, standard gravitational parameter
    ###        => G is the gravitational constant
    ###        => M is the mass of the more massive body (CLOCK, here)
    ###
    ### Solve for a and speed, assuming circular orbit and mu = 1.0:
    ###
    ###   a = (T / (2 PI))**(2/3)
    ###   v = (2 PI a) / T
    ####    = ((2 PI) / T) * (T / (2 PI))**(2/3)
    ####    = ((2 PI) / T)**(1/3)
    ####    = a**(-1/2)

    ### "Period" for minute hand is one hour
    aMinute = (sph/twopi) ** (2.0/3.0)
    spdMinute = aMinute ** (-0.5)

    ### "Period" for hour hand is half a day
    aHour = ((spd/2)/twopi) ** (2.0/3.0)
    spdHour = aHour ** (-0.5)
 
    ### Open new SP-Kernel (SPK; ephemerides of MINUTE and HOUR wrt CLOCK)
    handle = sp.spkopn(fnClockSpk, 'clock_simulation', 0)

    ### Equate the direction from the center of the clock towards the 12
    ### with the +X axis in the J2000 frame

    ### At the top of every hour, the minute hand is at [+aMinute,0,0],
    ### with a speed of spdMinute along +Y, i.e. velocty = [0,+spdMinute,0]

    ### At midnight, the hour hand is at [+aHour,0,0],
    ### with a speed of spdHour along +Y, i.e. velocty = [0,+spdHour,0]

    ### Store such states at (ET0 - one day) and (ET0 + two days)
    sp.spkw05( handle
             , iMinute, iClock, 'J2000'
             , -spd, 2*spd
             , 'minute_orbit'
             , 1.0, 2
             , [[aMinute,0.,0.,  0.,spdMinute,0.]
               ,[aMinute,0.,0.,  0.,spdMinute,0.]
               ]
             , [-spd,2*spd]
             )

    sp.spkw05( handle
             , iHour, iClock, 'J2000'
             , -spd, 2*spd
             , 'hour_orbit'
             , 1.0, 2
             , [[aHour,0.,0.,  0.,spdHour,0.]
               ,[aHour,0.,0.,  0.,spdHour,0.]
               ]
             , [-spd,2*spd]
             )

    ### Close SPK
    sp.spkcls(handle)

    ### End of clock SPK creation
    ####################################################################

  ### FURNSH the clock SPK
  sp.furnsh(fnClockSpk)

  ######################################################################
  ### Test the clock SPK
  halfHour = sph / 2.0
  degPerHour = 2.0 * dpr * twopi / hpd

  xTolerance = lambda xDiff: xDiff < 1e-10
  xDiffTolerance = lambda x,xExpect: xTolerance(abs(x-xExpect))

  et,iPass = et0, 0
  while et < (et0+spd+1):

    ### Get the state of the MINUTE and HOUR bodies every half hour
    stMinute,lt = sp.spkezr(sMinute, et, 'J2000', 'NONE', sClock)
    stHour,lt = sp.spkezr(sHour, et, 'J2000', 'NONE', sClock)

    if (iPass % 2):
      ### On the half hour, the minute hand will be at RA=180, along [-1,0,0]
      assert xTolerance(sp.vnorm(sp.vsub([-1.0,0.,0.], sp.vhat(stMinute[:3]))))
    else:
      ### On the hour, the minute hand will be at RA=0, at [+1,0,0]
      assert xTolerance(sp.vnorm(sp.vsub([1.0,0.,0.], sp.vhat(stMinute[:3]))))
      vsepDeg = dpr * sp.vsep(stHour[:3],stMinute[:3])
      vsepExpectDeg = (degPerHour * (iPass>>1)) % 360
      ### On the hour, the houre hand will be (30 * iPass) degrees clockwise
      ### from +X, and therefor also the same from the minute hand
      assert xDiffTolerance(vsepDeg, vsepExpectDeg) or xDiffTolerance((360 - vsepDeg), vsepExpectDeg)

    ### Step to next half hour and pass
    et += halfHour
    iPass += 1

  ### End of test
  ######################################################################

  ### Find classic problem:  how many times in a day the hour and minute
  ### hands are aligned

  cnfine = sp.stypes.SPICEDOUBLE_CELL(2)
  result = sp.stypes.SPICEDOUBLE_CELL(200)

  ### Set confinement window to 24h at 10ms before two successive midnights
  etStart,etStop = et0 - 10e-3, et0 + spd - 10e-3
  sp.wninsd(etStart, etStop, cnfine)

  ### Find local minima
  sp.gfpa(sClock, sMinute, "NONE",  sHour, "LOCMIN", 1e-6
         , 0.0, spm, 6000, cnfine, result)

  ### Confirm that 22 minima were found
  assert 22 == sp.wncard(result)

  ### Optional logging
  if doDebug:
    print('Alignments between {} and {}:'.format(sp.etcal(etStart+0.0005,99),sp.etcal(etStop+0.0005,99)))
    for iWin in range(sp.wncard(result)):
      left,right = sp.wnfetd(result,iWin)
      print('  {}:  {}'.format(iWin,sp.etcal(left+0.0005,99)))

  ### Repeat with extra time
  cnfine = sp.stypes.SPICEDOUBLE_CELL(2)
  etStop = et0 + spd + 10e-3
  sp.wninsd(etStart, etStop, cnfine)

  ### Find local minima
  sp.gfpa(sClock, sMinute, "NONE",  sHour, "LOCMIN", 1e-6
         , 0.0, spm, 6000, cnfine, result)

  ### Confirm that 23 minima were found
  assert 23 == sp.wncard(result)

  ### Optional logging
  if doDebug:
    print('Alignments between {} and {}:'.format(sp.etcal(etStart+0.0005,99),sp.etcal(etStop+0.0005,99)))
    for iWin in range(sp.wncard(result)):
      left,right = sp.wnfetd(result,iWin)
      print('  {}:  {}'.format(iWin,sp.etcal(left+0.0005,99)))

  ### Now search for events with a phase angle of 30deg
  cnfine = sp.stypes.SPICEDOUBLE_CELL(2)
  sp.wninsd(etStart, etStop, cnfine)

  ### Find local minima
  sp.gfpa(sClock, sMinute, "NONE",  sHour, "=", 30/dpr
         , 0.0, spm, 6000, cnfine, result)

  ### Confirm that 23 minima were found
  assert 44 == sp.wncard(result)

  ### Optional logging
  if doDebug:
    print('30deg events between {} and {}:'.format(sp.etcal(etStart+0.0005,99),sp.etcal(etStop+0.0005,99)))
    for iWin in range(sp.wncard(result)):
      left,right = sp.wnfetd(result,iWin)
      print('  {}:  {}'.format(iWin,sp.etcal(left+0.0005,99)))
