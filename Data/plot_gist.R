


par(mar = 0*c(5.1,4.1,4.1,2.1))

mybreaks= seq(min(dat_full$worldclim_temp_mid)-0.1, max(dat_full$worldclim_temp_mid)+0.1, length.out=100)
mycols =   snapalette( "Ipanema", length(mybreaks)-1, type = "continuous" )
mm = match(shp1$adm0_adm1,dat_full$adm0_adm1)
vcols = findInterval(dat_full$worldclim_temp_mid,mybreaks)



plot(shp0, xlim=c(-15,45),ylim=c(-20,35))
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey70",add=T) 
plot(shp1[!is.na(mm),],col=mycols[vcols], xlim=c(-15,45),ylim=c(-20,35) , lty=0, add=T)


image.plot(legend.only=TRUE,breaks=mybreaks,col=mycols,zlim=c(0,1), horizontal = TRUE, 
           legend.mar = 3.5)

########################

par(mar = 0*c(5.1,4.1,4.1,2.1))

mybreaks= seq(min(dat_full$worldclim_temp_range)-0.1, max(dat_full$worldclim_temp_range)+0.1, length.out=100)
mycols =  rev( snapalette( "Bouquet", length(mybreaks)-1, type = "continuous" ))
mm = match(shp1$adm0_adm1,dat_full$adm0_adm1)
vcols = findInterval(dat_full$worldclim_temp_range,mybreaks)



plot(shp0, xlim=c(-15,45),ylim=c(-20,35))
mm0 = match(shp0$ISO,c34) #
plot(shp0[is.na(mm0),],col="grey70",add=T) 
plot(shp1[!is.na(mm),],col=mycols[vcols], xlim=c(-15,45),ylim=c(-20,35) , lty=0, add=T)


image.plot(legend.only=TRUE,breaks=mybreaks,col=mycols,zlim=c(0,1), horizontal = TRUE, 
           legend.mar = 3.5)