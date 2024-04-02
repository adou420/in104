CC = gcc
CFLAGS = −Wall −Wextra −s t d=c99 −g
LDFLAGS = −lSDL2 −lm −lSDL2 image
SRC = SDL example . c
EXECUTABLE = output

all: $(EXECUTABLE)

$(EXECUTABLE): $(SRC)
	$(CC) $(CFLAGS) −o $@ $ˆ $(LDFLAGS)

run: $(EXECUTABLE)	
	./$(EXECUTABLE)

debug: $(EXECUTABLE)
	gdb ./$(EXECUTABLE)

clean:
	rm −f $(EXECUTABLE)