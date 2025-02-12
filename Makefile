NAME := bin/output

SRC_DIR := src
INC_DIR := include
OBJ_DIR := obj
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

CXX := g++
CXXFLAGS := -fdiagnostics-color=always -I$(INC_DIR) `root-config --cflags` -O3 -funroll-loops -Wall -Wextra -Wpedantic -Wshadow
LDFLAGS := `root-config --libs` -lvdt -flto=auto -Wl,--as-needed

$(NAME): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: debug
debug: CXXFLAGS := -fdiagnostics-color=always -I$(INC_DIR) `root-config --cflags` -O0 -Wall -Wextra -Wpedantic -Wshadow -g3
debug: LDFLAGS := `root-config --libs` -lvdt -Wl,--as-needed
debug: $(NAME)


.PHONY: clean
clean:
	rm -fv $(OBJS) $(NAME)

.PHONY: ASan
ASan: CXXFLAGS += -fsanitize=address -fno-omit-frame-pointer -g
ASan: LDFLAGS += -fsanitize=address
ASan: $(NAME)

.PHONY: gprof
gprof: CXXFLAGS += -pg
gprof: LDFLAGS += -pg
gprof: $(NAME)
