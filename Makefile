CC      = gcc

# Default: native optimizations + LTO
# Use `make PORTABLE=1` for portable builds (no -march=native)
ifdef PORTABLE
CFLAGS  = -O3 -Wall -Wextra -std=c11 -pthread -flto
else
CFLAGS  = -O3 -Wall -Wextra -std=c11 -pthread -march=native -flto
endif

LIBS    = -lm -pthread

SRCDIR  = src
OBJDIR  = obj
BINDIR  = bin

# Source files (add new modules here)
SRCS    = $(wildcard $(SRCDIR)/*.c)
OBJS    = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))
DEPS    = $(OBJS:.o=.d)

TARGET  = $(BINDIR)/mdl-repeat

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS) | $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(OBJDIR) $(BINDIR)

-include $(DEPS)
