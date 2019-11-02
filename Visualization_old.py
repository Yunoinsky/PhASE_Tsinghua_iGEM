# %% Import packages needed.
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from cv2 import cv2 as cv
import os
from PIL import Image
# %% The (color, shape) tuple for every point denoting the corresponding protein type.
plotDict = {2: ('y', '.'), 3: ('y', '.'), 4: ('g', '.'), 5: ('g', '.')}
# %% Define the function for plotting.


def do_plot(dictPointSet: dict, bx, by, bz, savePath: str):
    plotDictKeys = list(plotDict.keys())
    dictTypeToIndex = {}
    plotSet = []
    index = 0
    for ptype in plotDictKeys:
        plotSet.append([[], [], []])
        dictTypeToIndex[ptype] = index
        index += 1
    for position in list(dictPointSet.keys()):
        ptype = dictPointSet[position]
        if (ptype in plotDictKeys):
            ptypeIndex = dictTypeToIndex[ptype]
            plotSet[ptypeIndex][0].append(position[0])
            plotSet[ptypeIndex][1].append(position[1])
            plotSet[ptypeIndex][2].append(position[2])
    fig = plt.figure()
    ax = Axes3D(fig)
    for ptypeIndex in range(len(plotSet)):
        ptype = plotDictKeys[ptypeIndex]
        ax.scatter(plotSet[ptypeIndex][0], plotSet[ptypeIndex][1], plotSet[ptypeIndex]
                   [2], c=plotDict[ptype][0], marker=plotDict[ptype][1], s=5)
        ax.set_xlim(left=0, right=bx)
        ax.set_ylim(bottom=0, top=by)
        ax.set_zlim(bottom=0, top=bz)
    plt.savefig(savePath)
    plt.close()


# %% Convert log files to images.
imageFormat = '.png'
folderPath = input('Enter the folder path where the log files are located: ')
folderPath = folderPath.strip('/').strip('\\')+'/'
parameterFile = open(folderPath+'parameters.h', 'r')
parameterLines = parameterFile.readlines()
for parameterLine in parameterLines:
    if 'const int b[3]={' in parameterLine:
        break
parameterLine = parameterLine.strip(' ').strip(
    'const int b[3]={').strip('};//The size of the particles\' range of motion.\n')
parameterLine = parameterLine.split(',')
parameterLine = tuple(map(int, parameterLine))
fileNumber = 0
while True:
    fileNumber += 1
    fileName = folderPath + str(fileNumber).zfill(7) + '.csv'
    if os.path.exists(folderPath + 'fig' + str(fileNumber).zfill(7) + imageFormat):
        continue
    try:
        fileObject = open(fileName, 'r')
    except FileNotFoundError:
        fileNumber -= 1
        if input('File not found. Wait (w) or quit (q): ') in ['q', 'Q']:
            break
        else:
            continue
    fileLines = fileObject.readlines()
    del fileLines[0]
    del fileLines[-1]
    dictPlot = {}
    for fileLine in fileLines:
        fileLineSplited = fileLine.split(',')
        dictPlot[(int(fileLineSplited[2]), int(fileLineSplited[3]),
                  int(fileLineSplited[4]))] = int(fileLineSplited[1])
    do_plot(dictPlot, *parameterLine, folderPath +
            'fig' + str(fileNumber).zfill(7) + imageFormat)
    fileObject.close()
# %% Create the video file.
videoName = 'video.avi'
if fileNumber > 1:
    video_dir = folderPath+'fig_'+videoName
    img = Image.open(folderPath+'fig' +
                     '1'.zfill(7)+imageFormat)
    img_size = img.size
    fourcc = cv.VideoWriter_fourcc('M', 'J', 'P', 'G')
    videoWriter = cv.VideoWriter(video_dir, fourcc, 30, img_size)
    for i in range(fileNumber):
        im_name = os.path.join(
            folderPath, 'fig'+str(i+1).zfill(7)+imageFormat)
        frame = cv.imread(im_name)
        videoWriter.write(frame)
    print('Video \''+'fig_'+videoName+'\' created!\n')
    videoWriter.release()
