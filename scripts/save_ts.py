def save_ts(eta,xi,name):
    (tides_hr.ustar[:,eta,xi]*100*100).to_netcdf(f'ts_ustar_{name}.nc')
    tides_hr.Tstar[:,eta,xi].to_netcdf(f'ts_Tstar_{name}.nc')
    return

save_ts(592,348,'PI_front')
save_ts(867,419,'RonneNW')
save_ts(886,481,'RonneNE')
save_ts(823,499,'RonneHIR')
save_ts(758,389,'RonneEvISPocket')
save_ts(769,396,'RonneEvISOutfl')
save_ts(607,354,'PI_GL')
save_ts(601,350,'PI_trunk')


