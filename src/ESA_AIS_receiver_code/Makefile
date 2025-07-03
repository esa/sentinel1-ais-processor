CC = gcc
#CFLAGS = -O2 -Wall # Shows all warnings
LDFLAGS = -lm

SRC = AIS_receiver.c AIS_subroutines.c
OBJ = $(SRC:.c=.o)
TARGET = AIS_receiver

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)