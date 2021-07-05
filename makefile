NAMES := delaunay geometry delaunay_img

ifeq ($(OS), Windows_NT)
	RM := del

	SRCDIR := src
	OBJDIR := bin

	TARGET := main
else
	RM := rm

	SRCDIR := ./src
	OBJDIR := ./bin

	TARGET := ./main
endif

CC := gcc
CFLAGS := -Wall -Werror -Wextra -Wpedantic -std=c17 \
		-O2 -lm -Wno-sign-compare -lpthread -mavx2

INCDIR := $(SRCDIR)/include

OBJTARGS := $(addsuffix .o, $(NAMES) main)
OBJS := $(addprefix $(OBJDIR)/, $(OBJTARGS))
SRCS := $(addsuffix .c, $(addprefix $(SRCDIR)/, $(NAMES) main))
DEPS := $(addsuffix .h, $(addprefix $(INCDIR)/, $(NAMES)))

.PHONY: clobber run

build: $(OBJTARGS)
	$(CC) -o $(TARGET) $(OBJS) $(CFLAGS)

run: build
	$(TARGET)

dbg: build clobber
	gdb $(TARGET)

clobber:
	$(RM) $(OBJDIR)/*.o

$(OBJTARGS): %.o : $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $(OBJDIR)/$@ $< $(CFLAGS)
