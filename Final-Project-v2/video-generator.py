import re
import cv2
import os

import time
start_time = time.time()



image_folder = "C:\\Users\\samda\\Documents\\GitHub\\Heat-Transfer\\Final-Project-v2\\images"
output_folder = "C:\\Users\\samda\\Documents\\GitHub\\Heat-Transfer\\Final-Project-v2\\"
video_name = output_folder + "final-project.mp4"
fps = 5




###############################################################
#  Generate Video
###############################################################

def generate_video():
    def atoi(text):
        # https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        # https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]


    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    images.sort(key=natural_keys)
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    cv2.destroyAllWindows()
    video.release()

    print("\n--- Video Done ---")
    print("--- %.6f seconds ---" % (time.time() - start_time))




if __name__ == "__main__":
    generate_video()
