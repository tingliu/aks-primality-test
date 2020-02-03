CC = gcc
BUILD ?= perf
PARALLEL ?= true

CFLAGS = -Wall -Werror -Wextra
LDFLAGS = -lgmp

ifeq ($(BUILD),debug)
CFLAGS += -g
else
ifeq ($(BUILD),perf)
CFLAGS += -O3 -g -fno-omit-frame-pointer
else
ifeq ($(BUILD),release)
CFLAGS += -O3
endif
endif
endif

ifeq ($(PARALLEL),true)
CFLAGS += -fopenmp -DUSE_PARALLEL
LDFLAGS += -fopenmp
endif

TARGET = aks
SOURCES = aks.c
OBJECTS = $(SOURCES:%.c=%.o)

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJECTS)

test: $(TARGET)
	bash -c 'time ./$(TARGET) selected.txt'

clean:
	rm -f $(TARGET) $(OBJECTS)
