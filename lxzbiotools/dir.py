import os
import shutil
from pathlib import Path


#path of imgr
# path = 'D:\\BaiduNetdiskDownload\\newim\\'

# #path of folder
# folderPath = 'D:\\BaiduNetdiskDownload\\folderSort\\'

# peopleNumber = 61
# #new 61 folder numbers as sort_folder_number[61]
# sort_folder_number = [x for x in range(0,peopleNumber)]

# makedir 61 folders
'''
demo功能说明：
在folderPath处新建60个文件夹，
图片存储在path处
给每个文件夹分配150张图片（将9000张图片平均分配到60个文件夹）

Tips:
1: os.path.join(path1,path2...)
this function is used to combine the path,it returns a path which is 'path1/path2...'

2: os.makedirs(path)
this function is used to make a directory(new folder) in the path param

3: shutil.move(oldPath,newPath)
this function is used to move file from param1 to param 2

4: os.path.exists(path)
this function is used to check the filePath(param1) whether exists
'''

def MoveFile(file_dir: Path,
             output_dir: Path,
             dir_nums: int):
    file_dir = os.path.abspath(file_dir)
    folderPath = os.path.abspath(output_dir)
    sort_folder_number = [x for x in range(0,dir_nums+1)]
    for number in sort_folder_number:
        new_folder_path = os.path.join(folderPath,'%s'%number) #new_folder_path is ‘folderPath\number'
        if not os.path.exists(new_folder_path):
            os.makedirs(new_folder_path)
            print("new a floder named "+str(number)+'at the path of '+ new_folder_path)

    #the taget file list
    file_list = os.listdir(file_dir)
    folderNumber = 1
    print('there are '+str(len(file_list))+' files at the path of '+file_dir)
    for i in range(0,len(file_list)):
        old_file_path = os.path.join(file_dir,str(i))
        if os.path.isdir(old_file_path):
            '''if the path is a folder,program will pass it'''
            print('file does not exist ,path=' + old_file_path+' it is a dir' )
            pass
        elif not os.path.exists(old_file_path):
            '''if the path does not exist,program will pass it'''
            print('file does not exist ,path='+old_file_path)
            pass
        else:
            '''define the number,it decides how many files each directory process'''
            number = int(len(file_list)/dir_nums)
            if(i%number ==0):
                folderNumber +=1
            new_file_path = os.path.join(folderPath,'%s'%(folderNumber))
    if not os.path.exists(new_file_path):
        print('not exist path:'+new_file_path)
        break
    shutil.move(old_file_path,new_file_path)
    print('success move file from '+ old_file_path +' to '+new_file_path)