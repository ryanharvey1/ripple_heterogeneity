import subprocess
import os
def call_dlc(basepath):
    # os.system("conda activate DEEPLABCUT")
    config_path = 'D:/github/DeepLabCut/ojr-heath_ryan-2022-02-22/config.ymal'
    analyze_videos = 'deeplabcut.analyze_videos('+config_path+', [''Y:/OJRproject/OJR40/day7/test_211105_113610/Basler acA640-750uc (22583500)_20211105_113622038.avi''], save_as_csv=True)'
    subprocess.run('activate DEEPLABCUT && python && import deeplabcut &&'+ analyze_videos + '&& source deactivate', shell=False)
    
    # subprocess.run('source activate environment-name && "enter command here" && source deactivate', shell=True)

    # import deeplabcut
    # video_file = r'Y:\OJRproject\OJR40\day7\test_211105_113610\Basler acA640-750uc (22583500)_20211105_113622038.avi'
    # deeplabcut.analyze_videos(config_path, [video_file], save_as_csv=True)
    print(basepath)

call_dlc('basepath')
# deeplabcut.analyze_videos(config_path, ['fullpath/analysis/project/videos/reachingvideo1.avi'], save_as_csv=True)
