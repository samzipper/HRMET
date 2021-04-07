HRMET <- function(datetime, daylightSavings=1, longitude, latitude, 
                  Tair, SWin, u, ea, pa, 
                  Zair=2, Zu=2,
                  LAI, h, Tsurf, 
                  albSoil, albVeg, emissSoil, emissVeg) {
  # Author: Sam Zipper
  #         University of Victoria
  #         samuelczipper@gmail.com
  #
  # Any publications using HRMET should cite:
  #   Zipper, S.C. & S.P. Loheide II (2014). Using evapotranspiration to
  #       assess drought sensitivity on a subfield scale with HRMET, a high
  #       resolution surface energy balance model. Agricultural & Forest
  #		Meteorology 197: 91-102. DOI: 10.1016/j.agrformet.2014.06.009
  #
  # This work is licensed under the Creative Commons Attribution-
  # NonCommercial-ShareAlike 4.0 International License. To view a
  # copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/
  # or send a letter to Creative Commons, 444 Castro Street, Suite 900    ,
  # Mountain View, California, 94041, USA.
  #
  # This script uses an energy balance approach to calculate
  # evapotranspiration for a single point at a single instant. It can be run
  # over a gridded model domain (such as a field) using gridded input data,
  # but no interactions between neighboring cells occur. It has been
  # validated for corn ET in south-central Wisconsin grown on sandy and silt
  # loam soils. Please see the paper mentioned above for a full model
  # description.
  #
  # Necessary inputs are:
  # (1) Site information:
  #      -datetime = Date & time in Date format
  #      -daylightSavings = daylight savings time? [1=yes (summer), 0=no]
  #      -longitude = Longitude [decimal degrees]
  #      -latitude = Latitude [decimal degrees]
  #      -albSoil = Albedo of soil [0-1]
  #      -albVeg = Albedo of vegetation [0-1]
  #      -emissSoil = Emissivity of soil [0-1]
  #      -emissVeg = Emissivity of vegetation [0-1]
  # (2) Meteorological data:
  #      -Tair = Air temperature, degrees celsius
  #      -SWin = Incoming shortwave radiation [W/m2]
  #      -u = Wind speed [m/s]
  #      -ea = Air vapor pressure [kPa] (can be calculated from Tair and
  #            relative humidity)
  #      -pa = Atmospheric pressure [kPa]
  #      -Zair = height of air temperature measurements [m]
  #      -Zu = height of wind speed measurements [m]
  # (3) Canopy structure data:
  #      -LAI = Leaf Area Index [m2/m2]
  #      -h = Canopy height [m]
  # (4) Canopy surface temperature:
  #      -Tsurf = Canopy surface temperature [deg C]
  #
  # Here is an example that can run this:
  #  Site: WIBU-1, July 30 2012
  #ET <- HRMET(datetime=ymd_hm("2012-07-30 11:00"), daylightSavings=1,
  #            longitude=89.38420, latitude=43.2966, 
  #            Tair=27.25, SWin=731, u=1.85, ea=2.2834, pa=97.883, 
  #            Zair=3, Zu=10,
  #            LAI=4.04, h=2.54, Tsurf=25.62, 
  #            albSoil=0.105, albVeg=0.2, emissSoil=0.945, emissVeg=0.94)
  #
  # Last major update: December 2014
  # Rewritten from Matlab to R: June 2017
  
  # load lubridate - needed for some datetime processes
  require(lubridate)
  
  # define trig functions using degrees as input
  sind <- function(x) sin(x*pi/180)
  cosd <- function(x) cos(x*pi/180)
  asind <- function(x) asin(x)*180/pi
  acosd <- function(x) acos(x)*180/pi
  
  # First, check if any of your inputs are NaNs - if so, return an error
  if (is.na(datetime+longitude+latitude+Tair+SWin+u+ea+pa+LAI+h+Tsurf+albSoil+albVeg+emissSoil+emissVeg)) stop("Missing input data")
  
  # Error check: is wind measurement height < plant height?
  if (Zu < h) stop("Zu < h; provide wind speed measurement from different height")
  if (Zair < h) stop("Zair < h; provide air temperature measurement from different height")
  
  ## Data processing, define constants
  # make sure nothing is less thna or equal to 0 (important if you are selecting
  # your input data from a probability distribution, seems to come up for
  # wind speed from time to time).
  if (u < 0.1)    u <- 0.1
  if (ea < 0.01)  ea <- 0.01
  if (LAI < 0.01) LAI <- 0.01
  if (h < 0.01)   h <- 0.01
  
  # Define constants, unit conversions from input data, and any other
  # calculations that don't change based on x,y,Tsurf,z
  T_K <- Tsurf+273.16          # [degK] - surface temperature
  Tair_K <- Tair+273.16        # [degK] - air temperature, in Kelvin
  stefBoltz <- 5.67e-8         # [W m-2 K-4] - Stefan-Boltzmann Constant
  cP <- 29.3                   # [J mol-1 degC-1] - specific heat of air, from Campbell & Normal Table A1
  mW <- 0.018                  # [kg mol-1] - molecular mass of water, from C&N table A1
  mA <- 0.029                  # [kg mol-1] - molecular mass of air, from C&N table A1
  rhoW <- 1000                 # [kg m-3] - density of water, from C&N table A2
  rhoA <- (44.6*pa*273.15)/(101.3*Tair_K) # [mol m-3] - molar density of air, C&N eq. 3.3
  esA <- 0.611                  # [kPa] - constant 'a' for es equations, C&N pg 41
  esB <- 17.502                 # [-] - constant 'b' for es equations, C&N pg 41
  esC <- 240.97                 # [degC] - constant 'c' for es equations, C&N pg 41
  k <- 0.4                      # [-] - von Karman's constant
  year <- year(datetime)        # Separate out year
  julianDay <- yday(datetime)   # get julian day (1=Jan 1, etc.)
  julianTime <- hour(datetime)+minute(datetime)/60  # get julian time, in hours (0-24)
  
  ## Calculate Energy Balance
  ## Calculate R
  
  # Calculations for sun position from Campbell & Norman, chapter 11
  f <- 279.575+0.9856*julianDay
  eqTime <- (1/3600)*(-104.7*sin(f*(pi/180)) + 
                        596.2*sin(2*f*(pi/180)) + 
                        4.3*sin(3*f*(pi/180)) -
                        12.7*sin(4*f*(pi/180)) - 
                        429.3*cos(f*(pi/180)) - 
                        2*cos(2*f*(pi/180)) + 
                        19.3*cos(3*f*(pi/180)))   # equation of time (C&N 11.4)
  while (longitude > 7.5) {
    longitude <- longitude - 15 # figure out # of degrees west of standard meridian (if negative, that means you're east!)
  }
  LC <- -longitude/15                               # calculate longitude correction
  solarNoon <- 12 + daylightSavings - LC - eqTime   # [hrs] - calculate solar noon time
  solarDec <- asind(0.39785*sind(278.97 + 0.9856*julianDay +
                                   1.9165*sind(356.6+0.9856*julianDay)))        # [deg] - solar declination angle (C&N eq. 11.2)
  zenith <- acosd(sind(latitude)*sind(solarDec) +
                    cosd(latitude)*cosd(solarDec)*cosd(15*(julianTime-solarNoon)))  # [deg] - solar zenith angle (C&N eq. 11.1)
  
  # Calculate downwelling LW radiation, following approach of Crawford &
  # Duchon (1999). Includes cloudiness effects.
  solConst <- 1361.5    # [W m-2] - average annual solar constant. source = wikipedia
  m <- 35*cosd(zenith)*(1224*cosd(zenith)*cosd(zenith)+1)^(-0.5)  # calculate airmass number
  
  # look up table for G_Tw (delta) in Tw calculation, from Smith (1966) Table 1
  # Using only summer values here, but table also has spring/winter/fall
  #    (potential future improvement- select values based on julianDay)
  if (latitude < 10){
    G_Tw <- 2.80
  } else if (latitude < 20){
    G_Tw <- 2.70
  } else if (latitude < 30){
    G_Tw <- 2.98
  } else if (latitude < 40){
    G_Tw <- 2.92
  } else if (latitude < 50){
    G_Tw <- 2.77
  } else if (latitude < 60){
    G_Tw <- 2.67
  } else if (latitude < 70){
    G_Tw <- 2.61
  } else if (latitude < 80){
    G_Tw <- 2.24
  } else if (latitude < 90){
    G_Tw <- 1.94
  } else {
    stop("Latitude > 90 degrees")
  }
  
  Tdew <- esC*log(ea/esA)/(esB - log(ea/esA))        # [degC] - dewpoint temperature, from C&N eq. 3.14
  atmWater <- exp(0.1133 - log(G_Tw+1) + 0.0393*Tdew)  # atmospheric precipitable water, from Crawford & Duchon
  
  TrTpg <- 1.021 - 0.084*(m*(0.000949*pa*10+0.051))^0.5  # C&D Eq. 8 - corrections for Rayleigh scattering, absorption by permanent gases
  Tw <- 1 - 0.077*(atmWater*m)^0.3                       # C&D Eq. 9 - correction for absorption by water vapor
  Ta <- 0.935^m                                          # C&D Eq. 10 - correction for scattering by aerosols
  
  # Calculate LWin based on cloudiness
  Rso <- solConst*cosd(zenith)*TrTpg*Tw*Ta               # [W m-2] - clear sky shortwave irradiance
  clf <- max(c(0,min(c((1-SWin/Rso), 1))))                 # [-] - cloudiness fraction, from 0-1. Crawford & Duchon (1999)
  emissSky <- (clf + (1-clf)*(1.22+0.06*sin((month(datetime)+2)*pi/6))*(ea*10/Tair_K)^(1/7))    # [-] - emissivity based on cloud fraction, from Crawford & Duchon Eq. 20
  LWin <- emissSky*stefBoltz*(Tair_K^4)                  # [W m-2] - total incoming absorbed longwave radiation from atmosphere
  
  # Calculate LWout based on separate vegetation & soil components
  # (two-source model)
  fc <- 1-exp(-0.5*LAI)                                            # Norman et al. (1995) Eq. 3 - fractional plant cover based on LAI
  LWoutVeg <- emissVeg*stefBoltz*(T_K^4)*(1-exp(0.9*log(1-fc)))  # [W m-2] - outgoing LW from vegetation
  LWoutSoil <- emissSoil*stefBoltz*(T_K^4)*exp(0.9*log(1-fc))    # [W m-2] - outgoing LW from soil
  LWout <- LWoutVeg + LWoutSoil                                    # [W m-2] - total outgoing longwave radiation
  
  # Calculate total SWout as the sum of vegetation & soil components.
  # (two-source model)
  SWoutVeg <- SWin*(1-exp(0.9*log(1-fc)))*albVeg  # [W m-2] - outgoing shortwave radiation from vegetation canopy, based on amount of SW radiation reaching ground (Norman et al. (1995) Eq. 13))
  SWoutSoil <- SWin*exp(0.9*log(1-fc))*albSoil    # [W m-2] - outgoing shortwave radiation from soil surface, based on amount of SW radiation reaching ground (Norman et al. (1995) Eq. 13))
  SWout <- SWoutVeg + SWoutSoil                     # [W m-2] - outgoing SW radiation as sum of soil and vegetation components
  
  # Calculate net radiation budget (R)
  R <- SWin-SWout+LWin-LWout     # [W m-2] - net radiation at surface
  
  ## Calculate Ground heat flux (G) based on amount of radiation reaching the ground
  G <- 0.35*R*exp(0.9*log(1-fc))     # Eq. 13 - calculate G as 35% of R reaching soil (Norman et al., 1995)
  
  ## Iterative H calculation
  
  #Raupach (1994) z0m, d values as function of LAI, h
  Cw <- 2     # empirical coefficient
  Cr <- 0.3   # empirical coefficient
  Cs <- 0.003 # empirical coefficient
  Cd1 <- 7.5  # empirical coefficient
  uMax <- 0.3 # empirical coefficient
  subRough <- log(Cw)-1+(1/Cw) # roughness-sublayer influence function
  
  u_uh <- min(c(uMax, ((Cs+Cr*(LAI/2))^0.5)))       # ratio of u*/uh
  d <- h*(1-(1-exp(-sqrt(Cd1*LAI)))/sqrt(Cd1*LAI))  # [m] - zero-plane displacement height
  z0m <- h*(1-d/h)*exp(-k*(1/u_uh)-subRough)  # [m] - roughness length for momentum transfer
  kB1 <- 2.3          # kB^-1 factor from Bastiaansen SEBAL paper to convert from z0m to z0h kB1=2.3 means z0h = 0.1*z0m, which corresponds to C&N empirical equation
  z0h <- z0m/exp(kB1) # [m] - roughness length for heat transfer
  
  ############## Begin Iteration - positive stability ###################
  # Iterative solution to H, rH, etc. starting from very low values of H and positive stability
  zeta <- 0.5       # initial stability factor for diabatic correction (zeta from C&N sec 7.5) - >0 when surface cooler than air
  Hiter <- 0.5      # Hiter is placeholder for H during iterative process
  changePerc <- 0.5 # arbitrary starting value
  i <- 0            # starting i for iterations
  while (abs(changePerc)>0.001) {  # set convergence criteria here, in percent
    i <- i+1    # advance iteration number by one
    Hiter_old <- Hiter  # store value from previous iteration
    
    # calculate diabatic correction factors based on zeta. from C&N
    if (zeta > 0){
      # stable flow, C&N equation 7.27
      diabM <- 6*log(1+zeta)   # diabatic correction for momentum transfer
      diabH <- diabM  # diabatic correction for heat transfer
    } else {
      # unstable flow, C&N equation 7.26
      diabH <- -2*log((1+(1-16*zeta)^0.5)/2)
      diabM <- 0.6*diabH
    }
    
    #  calculate u*, gHa based on diabatic correction factors. from C&N
    uStar <- u*k/(log((Zu-d)/z0m)+diabM) # [m] - friction velocity C&N eq. 7.24
    rHa <- 1/((k^2)*u*rhoA/(1.4*((log((Zu-d)/z0m)+diabM)*(log((Zair-d)/z0h)+diabH)))) # [m2 s mol-1]
    
    # by including rExcess, we can ignore the difference between Tsurf and Taero
    rExcess <- mA*log(z0m/z0h)/(rhoA*k*uStar)   # [m2 s mol-1] - excess resistance, from Norman & Becker (1995)
    
    rHtot <- rHa+rExcess # [m2 s mol-1] - total resistance
    
    # Calculate H and zeta for next iteration
    Hiter <- cP*(Tsurf-Tair)/rHtot                        # sensible heat flux
    zeta <- -k*9.8*Zu*Hiter/(rhoA*cP*Tair_K*(uStar^3)) # updated zeta value
    
    # calculate the percent change between iterations
    if (i>2) {
      changePerc <- (Hiter-Hiter_old)/Hiter_old
    } else {
      changePerc <- 0.5  # Arbitrary initial value for changePerc on first iteration
    }
    
    if (i == 10000) stop("10000 iterations, will not converge")  # non-convergence scenario. occasionally happens when u is very very very low (or negative).
  }
  H_low <- Hiter   # [W m-2] - converged H for low starting values
  
  if (H_low == 0) H_low <- 0.02   # set to 0.02 so you don't divide by 0 later on
  
  rm(Hiter)
  rm(changePerc)
  rm(zeta)
  rm(i)
  
  
  ############## Begin Iteration - negative stability ###################
  # Repeat iterative solution to H & rH, starting from very high values of H and negative stability
  zeta <- -0.5      # initial stability factor for diabatic correction (zeta from C&N sec 7.5) - >0 when surface cooler than air
  Hiter <- -0.5     # Hiter is placeholder for H during iterative process
  changePerc <- 0.5 # arbitrary starting value
  i <- 0            # starting i for iterations
  while (abs(changePerc)>0.001) {  # set convergence criteria here, in percent
    i <- i+1    # advance iteration number by one
    Hiter_old <- Hiter  # store value from previous iteration
    
    # calculate diabatic correction factors based on zeta. from C&N
    if (zeta > 0){
      # stable flow, C&N equation 7.27
      diabM <- 6*log(1+zeta)   # diabatic correction for momentum transfer
      diabH <- diabM  # diabatic correction for heat transfer
    } else {
      # unstable flow, C&N equation 7.26
      diabH <- -2*log((1+(1-16*zeta)^0.5)/2)
      diabM <- 0.6*diabH
    }
    
    #  calculate u*, gHa based on diabatic correction factors. from C&N
    uStar <- u*k/(log((Zu-d)/z0m)+diabM) # [m] - friction velocity C&N eq. 7.24
    rHa <- 1/((k^2)*u*rhoA/(1.4*((log((Zu-d)/z0m)+diabM)*(log((Zair-d)/z0h)+diabH)))) # [m2 s mol-1]
    
    # by including rExcess, we can ignore the difference between Tsurf and Taero
    rExcess <- mA*log(z0m/z0h)/(rhoA*k*uStar)   # [m2 s mol-1] - excess resistance, from Norman & Becker (1995)
    
    rHtot <- rHa+rExcess # [m2 s mol-1] - total resistance
    
    # Calculate H and zeta for next iteration
    Hiter <- cP*(Tsurf-Tair)/rHtot                        # sensible heat flux
    zeta <- -k*9.8*Zu*Hiter/(rhoA*cP*Tair_K*(uStar^3)) # updated zeta value
    
    # calculate the percent change between iterations
    if (i>2) {
      changePerc <- (Hiter-Hiter_old)/Hiter_old
    } else {
      changePerc <- 0.5  # Arbitrary initial value for changePerc on first iteration
    }
    
    if (i == 10000) stop("10000 iterations, will not converge")  # non-convergence scenario. occasionally happens when u is very very very low (or negative).
  }
  H_high <- Hiter   # [W m-2] - converged H for low starting values
  
  if (H_high == 0) H_high <- 0.02   # set to 0.02 so you don't divide by 0 later on
  
  rm(Hiter)
  rm(changePerc)
  rm(zeta)
  rm(i)
  
  if (abs((H_low-H_high))/H_low <= 0.01){
    # if H_low and H_high are within 1# of each other, you've converged on a universal solution!
    H <- (H_low+H_high)/2               # [W m-2] - take the mean of your two H calculations for your sensible heat flux
  } else {  # if both H_low and H_high converge, but the numbers aren't close to each other, something weird happening - probably imaginary numbers in your solution for some reason (negative value somewhere)
    stop("H_low and H_high are too far apart!")
  }
  
  
  ## Calculate ET rate!
  LE <- R-G-H     # [W m-2] - latent heat flux rate as residual of energy budget
  lamda <- (2.495-(2.36e-3)*Tsurf)*mW*1000000    # [J mol-1] - latent heat of vaporization of water, dependent on temp [degC]. Formula B-21 from Dingman 'Physical Hydrology'
  ET <- (mW/rhoW)*(60*60)*1000*LE/lamda # [mm hr-1] - evaporation rate
  ET <- max(c(ET, 0))     # set equal to 0 if it calculates something below 0 (seems to happen occasionally in areas with very low veg cover)
  
  return(ET)
}