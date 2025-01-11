NAME := test

SRC_DIR := src
INC_DIR := include
OBJ_DIR := obj
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

CXX := g++
CXXFLAGS := -I$(INC_DIR) `root-config --cflags` -O3 -flto
LDFLAGS := `root-config --libs` -lvdt

$(NAME): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: debug
debug: CXXFLAGS += -g
debug: $(NAME)


.PHONY: clean
clean:
	rm -fv $(OBJS) $(NAME)