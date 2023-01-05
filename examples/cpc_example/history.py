import netCDF4 as nc


data=nc.Dataset('output/history.nc')

print('data', data['mass'][...])

##for i in range(1001)
#species_prop=['x', 'y', 'z', 'vx','vy', 'vz', 'charge', 'mass']
#for prop in species_prop:
#    print(prop, data[prop][...][:,-1])
##for var in data.variables:
##    print('var',var)
