# define our processor specific flags
export PLATFORM_FLAGS = -arch x86_64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
export LDFLAGS = -Wl -arch x86_64
