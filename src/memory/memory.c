#ifdef _WIN32

#include <windows.h>
#include <psapi.h>

void get_windows_memory(int *memory_kb)
{
    PROCESS_MEMORY_COUNTERS memInfo;
    GetProcessMemoryInfo(GetCurrentProcess(), &memInfo, sizeof(memInfo));
    *memory_kb = memInfo.WorkingSetSize / 1024; // Convert to KB
}

#else

#include <stdio.h>

void get_linux_memory(int *mem_kb)
{
    FILE *file;
    char line[256];
    int found = 0;

    // Initialize to -1 in case of error
    *mem_kb = -1;

    // Open /proc/self/status
    file = fopen("/proc/self/status", "r");
    if (file == NULL)
    {
        return;
    }

    // Read line by line
    while (fgets(line, sizeof(line), file))
    {
        // Look for the VmRSS line (resident memory)
        if (strncmp(line, "VmRSS:", 6) == 0)
        {
            // Extract the value in kB
            if (sscanf(line, "VmRSS: %d kB", mem_kb) == 1)
            {
                found = 1;
                break;
            }
        }
    }

    fclose(file);

    // If VmRSS not found, return -1
    if (!found)
    {
        *mem_kb = -1;
    }
}
#endif