CC = gcc
CFLAGS = -Wall -Werror -fopenmp

SRC_DIR = code_examples
INCLUDE_DIR = include
OUT_DIR = bin

INCLUDE_FILES := $(wildcard $(INCLUDE_DIR)/*.c)

ifeq ($(OS),Windows_NT)
	TARGET = console_app.exe
	RM = del /F /Q
else
	TARGET = my_program
	RM = rm -f
endif

all: $(OUT_DIR)/$(TARGET)

$(OUT_DIR)/$(TARGET):  $(SRC_DIR)/console_app.c $(INCLUDE_FILES) | $(OUT_DIR)
	$(CC) $(CFLAGS) -o $@ $^

$(OUT_DIR):
	mkdir $(OUT_DIR)

clean:
	$(RM) $(OUT_DIR)/$(TARGET)