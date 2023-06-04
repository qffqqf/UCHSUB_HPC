def merg_sort(arr):
    num = len(arr)
    num_d = len(arr)
    list_of_arr = [arr]
    count = 1
    while count > 0:
        list_of_arr_ = []
        count = 0
        for arr_ in list_of_arr:
            arr1, arr2, flag = devide(arr_)
            list_of_arr_.append([arr1, arr2])
            count = count + flag
        list_of_arr = list_of_arr_.copy()
    while len(list_of_arr) > 1:
        list_of_arr_ = []
        for pairs in list_of_arr:
            arr1, arr2 = pairs[0:2]
            list_of_arr_.append(sort(arr1, arr2))
        list_of_arr = list_of_arr_.copy()
    return list_of_arr


def devide(arr):
    num = len(arr)
    if num > 1:
        arr1 = arr[0:int(num/2)]
        arr2 = arr[int(num/2):-1]
        flag = 1
    else:
        flag = 0
    return arr1, arr2, flag

def sort(arr1, arr2):
    num1 = len(arr1)
    num2 = len(arr2)
    arr = []
    i = 0
    j = 0
    while i < num1 or j < num2:
        if arr1[i] < arr2[j]:
            arr.append(arr1[i])
            i +=1
        else:
            arr.append(arr2[j])
            j +=1
    return arr



arr = [3,8,1,6,0,2]
print(merg_sort(arr))