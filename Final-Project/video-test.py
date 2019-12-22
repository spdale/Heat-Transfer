import cv2
import os

image_folder = "C:\\Users\\samda\\Documents\\GitHub\\Heat-Transfer\\Final-Project\\test-images"
output_folder = "C:\\Users\\samda\\Documents\\GitHub\\Heat-Transfer\\Final-Project\\"
video_name = output_folder + "test.mp4"
fps = 5

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()