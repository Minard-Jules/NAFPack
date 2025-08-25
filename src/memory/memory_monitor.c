#include <stdio.h>
#include <string.h>

#ifdef _WIN32

#include <windows.h>
#include <psapi.h>

void get_memory(int *memory_kb) {
    // Initialize to -1 in case of error
    *memory_kb = -1;

    PROCESS_MEMORY_COUNTERS memInfo;
    GetProcessMemoryInfo(GetCurrentProcess(), &memInfo, sizeof(memInfo));
    *memory_kb = memInfo.WorkingSetSize / 1024;  // Convert to KB
}
#elif __linux__

void get_memory(int *memory_kb) {
    FILE *file;
    char line[256];
    int found = 0;

    // Initialize to -1 in case of error
    *memory_kb = -1;

    // Open /proc/self/status
    file = fopen("/proc/self/status", "r");
    if (file == NULL) {
        return;
    }

    // Read line by line
    while (fgets(line, sizeof(line), file)) {
        // Look for the VmRSS line (resident memory)
        if (strncmp(line, "VmRSS:", 6) == 0) {
            // Extract the value in kB
            if (sscanf(line, "VmRSS: %d kB", memory_kb) == 1) {
                found = 1;
                break;
            }
        }
    }

    fclose(file);
}
#else

void get_memory(int *memory_kb) {
    // Initialize to -1 in case of error
    *memory_kb = -1;

    printf("Memory monitoring not implemented for this OS.\n");
}

#endif

int get_memory_usage() {
    int memory_kb = -1;
    get_memory(&memory_kb);
    return memory_kb;
}