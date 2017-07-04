__author__ = 'Elena Maria Daniela Hindinger'

import numpy as np
import cv2
import os

batch_folder = r'I:/Elena H/DanioVision Preprocessing/Raw Videos' # ENTER LOCATION OF FOLDER CONTAINING ALL THE VIDEOS YOU WANT TO PROCESS
results_directory = r'I:/Elena H/DanioVision Preprocessing'  # ENTER SAVING LOCATION
filelist = os.listdir(batch_folder)

for i in filelist:
    video_path = os.path.join(batch_folder, i)
    basepath, filename = os.path.split(video_path)
    folder_name = filename[:-4]
    print ('Processing file ' + folder_name)

    # create path for saving videos
    output_folder = os.path.join(results_directory, folder_name)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    cap = cv2.VideoCapture(video_path)
    count = 1

    while (cap.isOpened()):
       # Capture frame-by-frame
       ret, frame = cap.read()

       # Our operations on the frame come here
       gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
       # print 'Saving frame ' + str(count)
       cv2.imwrite(os.path.join(output_folder, ("%s_frame_%d.tiff" % (folder_name, count))), frame)    # save frame as TIFF file
       count += 1

    cap.release()
    cv2.destroyAllWindows()

print ('All done.')