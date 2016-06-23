config.starGalaxyColumn = ""
filters = ('K')
config.magColumnMap = dict([(f,f+'_med') for f in filters])
config.magErrorColumnMap = dict([(f, f+'_err') for f in filters])
config.indexFiles = [
    'index-16042700.fits',
    'index-16042701.fits',
    'index-16042702.fits',
    'index-16042703.fits',
    'index-16042704.fits',
    ]
