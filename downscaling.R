library(gts)
library(scam)

# ESM
global = read_gts("input/ipsl-cm5a-lr_historical_to_zs_monthly_195001_200512.nc4")

# change prime meridian
global2 = global
longitude(global2) = longitude(global2, "center")

# observations
obs = read_gts("input/avhrr-only-v2-humboldt-n.198109-201706.nc4")
obs$grid$mask = mask(obs$x)

# change coordinates to observations is another option
obs2 = obs
longitude(obs2) = longitude(obs2, "left")

# spatial subsetting (ESM)
sim = subset(global2, grid=obs)
sim = subset(global2, grid=obs, expand=5)
# temporal subseting
sim = window(sim, start=start(obs))
obs = window(obs, end=end(sim))

# first regriding: bilinear by default
sim1 = regrid(object = sim, grid = obs)
# second regridding: akima 
sim2 = regrid(object = sim, grid = obs, method="akima", control=list(link="log"))
# to be implemented: knn, IDW, krigging
sim3 = regrid(object = sim, grid = obs, method="akima", extrap=TRUE)

par(mfrow=c(2,2))
plot(sim)
plot(sim1)
plot(sim2)
plot(sim3)

par(mfrow=c(2,1))
plot(obs)
plot(sim3)

sim3 = sim3 - 273.15 # change to celsius

par(mfrow=c(2,1))
plot(obs)
plot(sim3)

# compute climatologies
clim_obs = climatology(obs)
clim_sim = climatology(sim3)

# delta bias correction
sim_corr = sim3 - clim_sim + clim_obs

par(mfrow=c(1,3))
zlim = range(range(obs, na.rm=TRUE), range(sim3, na.rm=TRUE))
plot(obs, zlim=zlim)
plot(sim3, zlim=zlim)
plot(sim_corr, zlim=zlim)

# exercise: delta correction for RCP 2.6

# quantile mapping
quant_sim = quantile(sim3)
quant_obs = quantile(obs)

dat_sim = melt(quant_sim)
dat_obs = melt(quant_obs)

dat = merge(dat_obs, dat_sim, all = TRUE)
# dat = dat[complete.cases(dat), ]
mod = gam(sst ~ s(to), data=dat)

mod$pred = predict(mod, newdata=dat)


# projection
future = read_gts("input/ipsl-cm5a-lr_rcp26_to_zs_monthly_200601_21001231.nc4")
longitude(future) = longitude(future, "center")
fsim = subset(future, grid=obs, expand=5)
fsim = regrid(object = fsim, grid = obs, method="akima", extrap=TRUE)
fsim = fsim - 273.15
fdat = melt(fsim)
fdat$sst = predict(mod, newdata=fdat)
fsim_bc = fsim
fsim_bc$x[] = fdat$sst

par(mfrow=c(1,2))
plot(fsim)
plot(fsim_bc)
