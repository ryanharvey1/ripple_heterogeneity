import loading

basepath = r'Z:\Data\OMLproject\OML11\day16'
# basepath = r'Z:\Data\Kenji\ec013.152_157'

nChannels, fs, fs_dat, shank_to_channel = loading.loadXML(basepath)

print(fs)

print('loading lfp')
data, timestep = loading.loadLFP(
                                basepath,
                                n_channels=nChannels,
                                channel=1,
                                frequency=fs,
                                )
print(data.shape)
print(timestep.shape)

print('loading dat')
data, timestep = loading.loadLFP(
                                basepath,
                                n_channels=nChannels,
                                channel=1,
                                frequency=fs_dat,
                                ext='dat'
                                )
print(data.shape)
print(timestep.shape)

