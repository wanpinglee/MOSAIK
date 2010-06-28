# define our processor specific flags
export PLATFORM_FLAGS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Dstat64=stat -Dfstat64=fstat
export LDFLAGS = -Wl
