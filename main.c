/**
 * @file main.c
 * @author MC2018
 * @date 2022-03-20
 * @brief This project is for sorting a list via multiple threads
 * The approach taken for multithreading slightly varies from what was mentioned in the project, for efficiency purposes.
 * Instead of the threads having to wait every 50ms for the main method to restart the threads with new information,
 * the threads simply choose the next largest piece in their own logic, continuing until all pieces have been used up.
 * When finished, the main method's 50ms sleep timer then joins the threads as they finish up their work.
 * The threads call the receiveNextLargestPiece() method to be given the next piece.
 * If two threads attempt to run this method at the same time, one will be stalled until the other is finished in order to avoid collision.
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/time.h>
#include <pthread.h>

#define READ_END 0
#define WRITE_END 1
#define INT_SIZE 4

// random seed
#define NONE 0
#define RANDOM_SEED -1

// alternate sort
#define QUICK_SORT 0
#define SHELL_SORT 1
#define INSERTION_SORT 2

// argument definitions
#define SIZE "-n"
#define THRESHOLD "-s"
#define ALTERNATE "-a"
#define SEED "-r"
#define MULTITHREAD "-m"
#define PIECES "-p"
#define MAXTHREADS "-t"
#define MEDIAN "-m3"

// arguments
int size = NONE;
int threshold = 10;
int alternateAlgorithm = SHELL_SORT;
int seed = NONE;
bool isMultiThreaded = false;
int amtOfPieces = 10;
int maxThreads = 4;
bool median = false;

// list and piece management
int* list;
int piecesWorkedOnCount = 0;
int piecesGenerated = 0;
int* pieceBounds;
int* pieceBoundsCopy;
bool* piecesWorkedOn;
bool currentlySelectingPiece = false;

/**
 * @brief Calculates the CPU time difference between two clock_ts
 * 
 * @param c1 start clock_t
 * @param c2 end clock_t
 * @return float (time difference in seconds)
 */
float cpuTimeDiff(clock_t c1, clock_t c2) {
    return ((double) ((c2 - c1) / 1000)) / 1000;
}

/**
 * @brief Calculates the time difference between two struct timevals
 * 
 * @param tv1 start timeval
 * @param tv2 end timeval
 * @return float (time difference in seconds)
 */
float timeDiff(struct timeval tv1, struct timeval tv2) {

    float seconds = tv2.tv_sec - tv1.tv_sec;
    float millis = (float)(((int)(tv2.tv_usec - tv1.tv_usec)) / 1000) / 1000;

    if (millis < 0) {
        millis *= -1;
    }

    return seconds + millis;
}

/**
 * @brief Sets up the random number generator seed
 * 
 */
void determineRandom() {
    if (seed == RANDOM_SEED) {
        srand(time(NULL) % 1000000);
    } else {
        srand(seed);
    }
}

/**
 * @brief Runs a shell sort on the list
 * 
 * @param low (low index)
 * @param highPosition (high index + 1)
 */
void shellsort(int low, int highPosition) {
    int p = 0, k = 1, temp;
    int n = highPosition + 1 - low;

    while (k <= n) {
        k *= 2;
    }

    k = (k / 2) - 1;

    do {
        for (int i = low; i < (highPosition - k); i++) {
            for (int j = i; j >= low; j -= k) {
                if (list[j] <= list[j + k]) {
                    break;
                }
                
                temp = list[j];
                list[j] = list[j + k];
                list[j + k] = temp;
            }
        }

        k /= 2;
    } while (k > 0);
}

/**
 * @brief Selects the middlemost element between a low, middle, and high index
 * 
 * @param low (lowest index to read from)
 * @param high (highest index to read from)
 */
void medianOfThree(int low, int high) {
    int middle = (high + low) / 2;
    int temp, median, l = list[low], m = list[middle], h = list[high];

    if ((l >= m && l <= h) || (l >= h && l <= m)) {
        median = low;
    } else if ((m >= l && m <= h) || (m >= h && m <= l)) {
        median = middle;
    } else {
        median = high;
    }

    temp = list[low];
    list[low] = list[median];
    list[median] = temp;
}

/**
 * @brief Runs an insertion sort on the list
 * 
 * @param low (low index)
 * @param highPosition (high index + 1)
 */
void insertionsort(int low, int highPosition) {
    int i = low, j, temp;
    
    while (i < highPosition) {
        j = i;

        while (j > low && list[j - 1] > list[j]) {
            temp = list[j];
            list[j] = list[j - 1];
            list[j - 1] = temp;
            j--;
        }

        i++;
    }
}

/**
 * @brief Partitions a piece of the list -- used for initialization of pieces
 * There's low enough overhead for pieces that it won't greatly impact performance, so it was put in a separate method
 * 
 * @param low (lowest index)
 * @param high (highest index)
 * @return int (index of partition)
 */
int partition(int low, int high) {
    int pivot = list[high];
    int i = low - 1, temp;
    
    if (median) {
        medianOfThree(low, high);
    }

    for (int j = low; j < high; j++) {
        if (list[j] <= pivot) {
            temp = list[++i];
            list[i] = list[j];
            list[j] = temp;
        }
    }

    temp = list[i + 1];
    list[i + 1] = list[high];
    list[high] = temp;
    return i + 1;
}

/**
 * @brief Runs the customized quicksort on the list
 * For each of the two sections (low to pivot and pivot to high), it checks for several cases:
 * 1. If the size of the section spanned is 1, do nothing
 * 2. If the size of the section spanned is 2, check to see if they need swapped
 * 3. If the size of the section spanned is greater than the threshold, do another quicksort
 * 4. Otherwise, choose between insertion sort and shell sort based on what is selected
 * 
 * @param low (low index)
 * @param high (high index)
 */
void quicksort(int low, int high) {
    if (low < high) {
        int highPivot = list[high];
        int pivot = low - 1, temp, spanSize;
    
        if (median) {
            medianOfThree(low, high);
        }

        for (int j = low; j < high; j++) {
            if (list[j] <= highPivot) {
                temp = list[++pivot];
                list[pivot] = list[j];
                list[j] = temp;
            }
        }

        temp = list[pivot + 1];
        list[pivot + 1] = list[high];
        list[high] = temp;
        pivot++;

        // quicksort logic
        spanSize = pivot - low;
        
        if (spanSize <= 1) {
        } else if (spanSize == 2) {
            if (list[low] > list[low + 1]) {
                int temp = list[low];
                list[low] = list[low + 1];
                list[low + 1] = temp;
            }
        } else if (spanSize > threshold) {
            quicksort(low, pivot - 1);
        } else if (alternateAlgorithm == SHELL_SORT) {
            shellsort(low, pivot);
        } else if (alternateAlgorithm == INSERTION_SORT) {
            insertionsort(low, pivot);
        }

        spanSize = high - pivot;

        if (spanSize == 1) {
        } else if (spanSize == 2) {
            if (list[pivot + 1] > list[high]) {
                int temp = list[pivot + 1];
                list[pivot + 1] = list[high];
                list[high] = temp;
            }
        } else if (spanSize > threshold) {
            quicksort(pivot + 1, high);
        } else if (alternateAlgorithm == SHELL_SORT) {
            shellsort(pivot + 1, high + 1);
        } else if (alternateAlgorithm == INSERTION_SORT) {
            insertionsort(pivot + 1, high + 1);
        }
    }
}

/**
 * @brief Generates and randomizes the list
 * 
 */
void generateRandomList() {
    struct timeval start, end;
    gettimeofday(&start, NULL);

    list = malloc(size * sizeof(int));

    gettimeofday(&end, NULL);
    printf("Array created in \t%f\n", timeDiff(start, end));
    gettimeofday(&start, NULL);

    for (int i = 0; i < size; i++) {
        list[i] = i;
    }

    gettimeofday(&end, NULL);
    printf("Array initialized in \t%f\n", timeDiff(start, end));
    gettimeofday(&start, NULL);

    for (int i = 0; i < size; i++) {
        int random = rand() % size;
        int temp = list[i];
        list[i] = list[random];
        list[random] = temp;
    }

    gettimeofday(&end, NULL);
    printf("Array randomized in \t%f\n", timeDiff(start, end));
}

/**
 * @brief Reads/assigns values for all arguments and checks to make sure arguments are valid
 * 
 * @param argc (size of argument list)
 * @param argv (list of arguments)
 * @return true if no failures arise
 * @return false if a problem with an argument occurs
 */
bool determineArguments(int argc, char* argv[]) {
    int i;
    bool errorFound = false;

    for (i = 0; i < argc; ) {
        if (i + 1 == argc) {
            i = argc;
            continue;
        }

        i++;

        if (strcmp(argv[i - 1], SIZE) == 0) {
            size = atoi(argv[i]);

            if (size <= 0) {
                errorFound = true;
                break;
            }
        } else if (strcmp(argv[i - 1], THRESHOLD) == 0) {
            threshold = atoi(argv[i]);

            if (threshold <= 0 && argv[i][0] != '0') {
                errorFound = true;
                break;
            }
        } else if (strcmp(argv[i - 1], ALTERNATE) == 0) {
            if (strcmp(argv[i], "s") == 0 || strcmp(argv[i], "S") == 0) {
                alternateAlgorithm = SHELL_SORT;
            } else if (strcmp(argv[i], "i") == 0 || strcmp(argv[i], "I") == 0) {
                alternateAlgorithm = INSERTION_SORT;
            } else {
                errorFound = true;
                break;
            }
        } else if (strcmp(argv[i - 1], SEED) == 0) {
            seed = atoi(argv[i]);

            if (strcmp(argv[i], "-1") == 0) {
                seed = RANDOM_SEED;
            } else if (seed <= 0 && argv[i][0] != '0') {
                errorFound = true;
                break;
            }
        } else if (strcmp(argv[i - 1], MULTITHREAD) == 0) {
            if (strcmp(argv[i], "n") == 0 || strcmp(argv[i], "N") == 0) {
                isMultiThreaded = false;
            } else if (strcmp(argv[i], "y") == 0 || strcmp(argv[i], "Y") == 0) {
                isMultiThreaded = true;
            } else {
                errorFound = true;
                break;
            }
        } else if (strcmp(argv[i - 1], PIECES) == 0) {
            amtOfPieces = atoi(argv[i]);

            if (amtOfPieces <= 0) {
                errorFound = true;
                break;
            }
        } else if (strcmp(argv[i - 1], MAXTHREADS) == 0) {
            maxThreads = atoi(argv[i]);

            if (maxThreads <= 0) {
                errorFound = true;
                break;
            }
        } else if (strcmp(argv[i - 1], MEDIAN) == 0) {
            if (strcmp(argv[i], "n") == 0 || strcmp(argv[i], "N") == 0) {
                median = false;
            } else if (strcmp(argv[i], "y") == 0 || strcmp(argv[i], "Y") == 0) {
                median = true;
            } else {
                errorFound = true;
                break;
            }
        }
    }

    if (errorFound) {
        printf("Improper syntax for %s: %s\n", argv[i - 1], argv[i]);
        return true;
    }

    if (size == NONE) {
        printf("No size (%s) was specified.\n", SIZE);
        return true;
    }

    if (!isMultiThreaded) {
        maxThreads = 1;
    }

    if (maxThreads > amtOfPieces) {
        maxThreads = amtOfPieces;
    }

    return false;
}

/**
 * @brief Checks to verify whether or not the algorithm sorted correctly
 * 
 * @return int (# of incorrectly ordered numbers)
 */
int isSorted() {
    int wrong = 0;

    for (int i = 1; i < size; i++) {
        if (list[i - 1] + 1 != list[i]) {
            wrong++;
        }
    }

    return wrong;
}

/**
 * @brief Finds the largest piece's starting index
 * 
 * @return int (starting index of the largest piece)
 */
int findLargestPieceStartingIndex() {
    int maxIndex = 0;

    for (int i = 2; i < piecesGenerated * 2; i += 2) {
        if (pieceBounds[i] == -1) {
            continue;
        }

        if (pieceBounds[i + 1] - pieceBounds[i] > pieceBounds[maxIndex + 1] - pieceBounds[maxIndex]) {
            maxIndex = i;
        }
    }

    return maxIndex;
}

/**
 * @brief Generates the pieces from the list
 * 
 */
void generatePieces() {
    pieceBounds = malloc(amtOfPieces * sizeof(int) * 2);
    pieceBoundsCopy = malloc(amtOfPieces * sizeof(int) * 2);
    piecesWorkedOn = malloc(amtOfPieces * sizeof(int));
    pieceBounds[0] = 0;
    pieceBounds[1] = size - 1;
    piecesGenerated++;
    piecesWorkedOnCount = 0;

    for (int i = 0; i < amtOfPieces; i++) {
        piecesWorkedOn[i] = false;
    }

    if (amtOfPieces == 1) {
        return;
    }

    while (piecesGenerated != amtOfPieces) {
        int largestStartingIndex = findLargestPieceStartingIndex();
        int partitionIndex = partition(pieceBounds[largestStartingIndex], pieceBounds[largestStartingIndex + 1]);

        // These two sets of if statements are to ensure safe sorting on small-sized sorts
        if (partitionIndex == pieceBounds[largestStartingIndex + 1]) {
            pieceBounds[piecesGenerated * 2] = partitionIndex;
        } else {
            pieceBounds[piecesGenerated * 2] = partitionIndex + 1;
        }

        pieceBounds[piecesGenerated * 2 + 1] = pieceBounds[largestStartingIndex + 1];

        if (partitionIndex == pieceBounds[largestStartingIndex]) {
            pieceBounds[largestStartingIndex + 1] = partitionIndex;
        } else {
            pieceBounds[largestStartingIndex + 1] = partitionIndex - 1;
        }

        piecesGenerated++;
    }

    for (int i = 0; i < amtOfPieces * 2; i++) {
        pieceBoundsCopy[i] = pieceBounds[i];
    }
}

/**
 * @brief Marks the largest non-worked on piece as worked on
 * 
 * @return int (starting index of the largest piece yet to be worked on)
 */
int markLargestPieceWorkedOn() {
    int maxStartingIndex = 0;

    for (int i = 2; i < piecesGenerated * 2; i += 2) {
        if (pieceBoundsCopy[i] == -1) {
            continue;
        }

        if (pieceBoundsCopy[i + 1] - pieceBoundsCopy[i] > pieceBoundsCopy[maxStartingIndex + 1] - pieceBoundsCopy[maxStartingIndex]) {
            maxStartingIndex = i;
        }
    }

    pieceBoundsCopy[maxStartingIndex] = -1;
    pieceBoundsCopy[maxStartingIndex + 1] = -1;
    return maxStartingIndex;
}

/**
 * @brief Finds the next largest piece to work on and marks it to not be read again
 * 
 * @return int (next largest piece starting index)
 */
int receiveNextLargestPiece() {
    int startingIndex = 0;

    // Used as a waiting system to make sure no two threads receive the same piece
    while (currentlySelectingPiece) {
        usleep(100);
    }

    currentlySelectingPiece = true;

    piecesWorkedOnCount++;
    startingIndex = markLargestPieceWorkedOn();

    currentlySelectingPiece = false;
    return startingIndex;
}

/**
 * @brief The method run as a thread to sort the list
 * 
 * @param param 
 * @return void* 
 */
void *sortingThread(void *param) {
    int num;

    while (piecesWorkedOnCount != amtOfPieces) {
        num = receiveNextLargestPiece();
        
        quicksort(pieceBounds[num], pieceBounds[num + 1]);
    }

    pthread_exit(0);
}

/**
 * @brief Main method which runs the program -- self-explanatory
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, char* argv[]) {
    bool* finishedThreads;
    bool test;
    struct timeval totalTimeStart, totalTimeEnd; gettimeofday(&totalTimeStart, NULL);
    struct timeval sortingTimeStart, sortingTimeEnd;
    clock_t cpuSortingTimeStart, cpuSortingTimeEnd;
    bool argumentErrorDetected = determineArguments(argc, argv);
    int finishedThreadCounter = 0;
    pthread_t tids[maxThreads];
    pthread_attr_t attrs[maxThreads];

    if (argumentErrorDetected) {
        printf("Terminating program.\n");
        return 0;
    }

    finishedThreads = malloc(sizeof(bool) * maxThreads);

    determineRandom();
    generateRandomList();
    gettimeofday(&sortingTimeStart, NULL);
    cpuSortingTimeStart = clock();
    generatePieces();

    for (int i = 0; i < maxThreads; i++) {
        pthread_create(&tids[i], NULL, sortingThread, NULL);
        finishedThreads[i] = false;
    }

    while (finishedThreadCounter != maxThreads) {
        for (int i = 0; i < maxThreads; i++) {
            int s;

            if (finishedThreads[i]) {
                continue;
            }

            s = pthread_tryjoin_np(tids[i], NULL);

            if (s == 0) {
                finishedThreads[i] = true;
                finishedThreadCounter++;
            }
        }

        usleep(50000);
    }

    cpuSortingTimeEnd = clock();
    gettimeofday(&sortingTimeEnd, NULL);
    gettimeofday(&totalTimeEnd, NULL);

    printf("Seconds spent sorting: Wall Clock:\t %f / ", timeDiff(sortingTimeStart, sortingTimeEnd));
    printf("CPU: %f\n", cpuTimeDiff(cpuSortingTimeStart, cpuSortingTimeEnd));
    printf("Total time: \t%f\n", timeDiff(totalTimeStart, totalTimeEnd));
    printf("# Wrong: %d\n", isSorted());
    return 0;
}