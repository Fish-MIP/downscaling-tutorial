install.packages(c("akima", "fields", "lubridate", "maps", "maptools", "methods",
                   "mgcv", "ncdf4", "sp", "remotes"))
remotes::install_github("cran/maptools")
remotes::install_github("roliveros-ramos/colorful")
remotes::install_github("roliveros-ramos/nctools")
remotes::install_github("roliveros-ramos/gts")

library(gts)

# ESM
global = read_gts("input/ipsl-cm5a-lr_historical_to_zs_monthly_195001_200512.nc4")

# change prime meridian
global2 = global
longitude(global2) = longitude(global2, "center")

global3 = global2
longitude(global3) = longitude(global3, "left")

                                                                                                                                                                                                                                                                                                                                                                  
# observations
obs = read_gts("input/avhrr-only-v2-humboldt-n.198109-201706.nc4")
# obs$grid$mask = mask(obs$x)
mask(obs) = mask(obs$x)
# other_mask = mask(obs)

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

dif = obs - sim3
col = colorful::divergencePalette(
  zlim=range(dif, na.rm=TRUE), symmetric = TRUE)
plot(dif, col=col)

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

range(obs, na.rm=TRUE)
range(sim3, na.rm=TRUE)
range(sim_corr, na.rm=TRUE)

mean(obs, na.rm=TRUE)
media_tiempo = mean(obs, by="time")
media_ori = mean(sim3, by="time")
media_sim = mean(sim_corr, by="time")
media_espacio = mean(obs, by="space")

par(mfrow=c(2,1))

plot(media_tiempo)
lines(media_sim, col="red", lwd=2)
lines(media_ori, col="blue", lwd=2)

image(media_espacio)

# exercise: delta correction for RCP 2.6

# leer escenarios
rcp26 = read_gts("input/ipsl-cm5a-lr_rcp26_to_zs_monthly_200601_21001231.nc4")
rcp85 = read_gts("input/ipsl-cm5a-lr_rcp85_to_zs_monthly_200601_21001231.nc4")

longitude(rcp26) = longitude(rcp26, "center")
longitude(rcp85) = longitude(rcp85, "center")
# subset
sim26 = subset(rcp26, grid=obs, expand=5)
sim85 = subset(rcp85, grid=obs, expand=5)

sim26 = regrid(object = sim26, grid = obs, method="akima", extrap=TRUE)
sim85 = regrid(object = sim85, grid = obs, method="akima", extrap=TRUE)

sim26 = sim26 - 273.15
sim85 = sim85 - 273.15

sim26 = sim26 - clim_sim + clim_obs
sim85 = sim85 - clim_sim + clim_obs

m26 = mean(sim26, by="time")
m85 = mean(sim85, by="time")

plot(media_sim, xlim=c(1980, 2100), ylim=c(20, 30), las=1)
lines(m26, col="blue")
lines(m85, col="red")

slice80.26 = window(sim26, start=c(2080,1), end=c(2089,12))
slice80.85 = window(sim85, start=c(2080,1), end=c(2089,12))

clim80.26 = climatology(slice80.26)
clim80.85 = climatology(slice80.85)

c80.26 = mean(clim80.26, by="time")
c80.85 = mean(clim80.85, by="time")

plot(mean(clim_obs, by="time"), ylim=c(20, 30), axes=FALSE)
lines(c80.26, col="blue")
lines(c80.85, col="red")
axis(2, las=1)
axis(1, at=(1:12 - 0.5)/12, labels=month.abb)
box()


# quantile mapping
probs = seq(0, 1, by=0.05)
quant_sim = quantile(sim3, probs=probs)
quant_obs = quantile(obs, probs=probs)

dat_sim = melt(quant_sim)
dat_obs = melt(quant_obs)

dat = merge(dat_obs, dat_sim, all = TRUE)
# dat = dat[complete.cases(dat), ]
library(mgcv)
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
zlim = range(range(fsim, na.rm=TRUE), range(fsim_bc, na.rm=TRUE))
plot(fsim, zlim=zlim)
plot(fsim_bc, zlim=zlim)

library(scam)

write_ncdf(fsim_bc, filename="output.nc")

