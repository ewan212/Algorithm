
'''
Binary search: 
- searching algo used in a sorted array 
- repeatedly divide the search interval in half
- O(logn)
'''




def iter_b_search(arr, x):
    '''
    iterative implementation of binary search
    :param: array to be searched from (will be sorted in the function)
    :x: search term
    '''
    arr.sort()
    low = 0
    high = len(arr) - 1
    while low <= high:
        mid = (low + high) // 2 

        if x == arr[mid]:
            return mid
            break
            
        elif x > arr[mid]:
            low = mid + 1
            
        else:
            high = mid - 1 
    return -1
            

def recursion_b_search(arr, x, low, high):
    '''
    recursive implementation of binary search
    '''
    arr.sort()
    if high >=0 :
        mid = (low + high) // 2 
        if x == arr[mid]:
            return mid
             
        elif x > arr[mid]:
            return recursion_b_search(arr, x, mid+1, high)
        else: 
            return recursion_b_search(arr, x, low, mid-1)
    else:
        return -1

arr = [2,4,10,9,8,45]
x = 2
recursion_b_search(arr, x, 0, len(arr) - 1 )