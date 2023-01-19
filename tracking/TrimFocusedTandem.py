# Video Cropping codes for Diacamma tandem run
# N. Mizumoto
# (initiate 7/2/2021, edit 10/11/2021, last check 6/18/2022)

# Note
# This code load original video and coordinate data for tandem run (created by FastTrack and then converted by TrackingConverter.R).
# Then will create a background to mask other area than focal tandem run.
# output will be in .mp4


# Loading Libraries
import cv2
import sys
import numpy as np
import glob
import csv
import pprint

## Loading data
HOME = 'the directly to the VideoPrepPython'
path = glob.glob(HOME + "videos/" + "*.mp4")
coordinate = glob.glob(HOME + "coordinate/" + "*.csv")
print(path)
print(coordinate)
filenums = list(range((len(path))))
# parameter for the area around the center of tandem
around = 20 #15

for i in filenums:
    
    #################################
    #### 1. Get file information ####
    #################################

    v = path[i]
    c = coordinate[i]
    print(v)
    print(c)

    ## read video
    video = cv2.VideoCapture(v)

    ## video information
    width = video.get(cv2.CAP_PROP_FRAME_WIDTH)
    height = video.get(cv2.CAP_PROP_FRAME_HEIGHT)
    count = video.get(cv2.CAP_PROP_FRAME_COUNT)
    fps = video.get(cv2.CAP_PROP_FPS)
    print("width:{}, height:{}, count:{}, fps:{}".format(width,height,count,fps))

    ## location information
    location = [[0 for i3 in range(3)] for i2 in range(10401)]
    with open(c) as f:
        reader = csv.reader(f)
        l = [row for row in reader]
        ln = np.asfarray(l, float)

    location = ln
    print(location)

    ############################
    #### 2. Make background ####
    ############################
    video = cv2.VideoCapture(v)
    has_next, i_frame = video.read()

    # background
    back_frame = i_frame.astype(np.float32)

    # converting
    num = 1
    while has_next == True:
        if num % 100 == 0:
            # convert image to float
            f_frame = i_frame.astype(np.float32)

            # diff calcuration
            diff_frame = cv2.absdiff(f_frame, back_frame)

            # update background
            cv2.accumulateWeighted(f_frame, back_frame, 0.025)

            # display the frame
            #cv2.imshow("WINDOW_BACK", back_frame.astype(np.uint8))
            #cv2.waitKey(1)
            #print(num)

        num = num + 1
        # next frame
        has_next, i_frame = video.read()

    cv2.imshow("WINDOW_BACK", back_frame.astype(np.uint8))
    cv2.waitKey(1)
    cv2.imwrite(HOME + "videos/" + "background.jpeg", back_frame)

    #########################
    #### 3. Video Create ####
    #########################
    fourcc = cv2.VideoWriter_fourcc('m','p','4', 'v')
    writer = cv2.VideoWriter(v.replace('.mp4', '_extract.mp4'), fourcc, fps, (int(width), int(height)))
    num = 0
    video.set(cv2.CAP_PROP_POS_FRAMES, 0) 

    video = cv2.VideoCapture(v)

    while(video.isOpened()):
        ret, frame = video.read()
        if num < location.shape[0]:
            x = int(location[num][1])
            y = int(location[num][2])

            #cv2.circle(frame, (x,y), 2, color=(0,0,255), thickness=-1)
            bg = cv2.imread(HOME + "videos/" + "background.jpeg")
            cv2.rectangle(bg, (x-around+1,y-around+1), (x+around-1, y+around-1), (0,0,0), -1, cv2.LINE_4)

            cv2.rectangle(frame, (1,1), (x-around, int(height)), (0,0,0), -1, cv2.LINE_4)
            cv2.rectangle(frame, (x+around,1), (int(width), int(height)), (0,0,0), -1, cv2.LINE_4)
            cv2.rectangle(frame, (1,1), (int(width), y-around), (0,0,0), -1, cv2.LINE_4)
            cv2.rectangle(frame, (1,y+around), (int(width), int(height)), (0,0,0), -1, cv2.LINE_4)
            
            frame=cv2.add(bg,frame)
            
            ## plot
            cv2.imshow('frame', frame)
            cv2.waitKey(1)

            writer.write(frame)
            #cv2.imwrite(HOME + "Plots/Visualize/p8-13/" + "picture{:0=5}".format(num)+".jpg", frame)
            #print("save picture{:0=5}".format(num)+".jpg")

        else:
            break

        num += 1

    cv2.destroyAllWindows()
    writer.release()

cv2.destroyAllWindows()