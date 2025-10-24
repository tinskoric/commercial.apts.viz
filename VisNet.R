source("Libraries.R")
library(igraph)
library(visNetwork)

# Formatting Note: data I format with underscores to split words;
# functions I format with capitalization to split words.

# data <- read_xls("data.xls", sheet = 1) %>% 
#  mutate(
#    from = ID,
#    to = ifelse(
#      Internal_Intersections == "NA",
#      ifelse(
#        External_Intersections == "NA",
#        NA,
#        External_Intersections
#      ),
#      ifelse(
#        External_Intersections == "NA",
#        Internal_Intersections,
#        paste0(Internal_Intersections, ", ", External_Intersections)
#      )
#    )
#  ) %>% separate_rows(to, sep = ",\\s*") %>% 
#  mutate(
#    type = ifelse(Internal_Intersections == to, "Internal", "External")
#  ) %>% mutate(
#    # below is to deal with the none cases in the network
#    to = ifelse(is.na(type), from, to),
#    type = ifelse(is.na(type), "None", type)
#  ) %>% select(from, to, type)
#
# data_igraph <- graph_from_data_frame(data, directed = TRUE)

################################################################################
# Function to format provided data to an igraph

Associations <- function(data) {
  data_igraph <- data %>% 
    mutate(
      from = ID,
      to = ifelse(
        Internal_Intersections == "NA",
        ifelse(
          External_Intersections == "NA",
          NA,
          External_Intersections
        ),
        ifelse(
          External_Intersections == "NA",
          Internal_Intersections,
          paste0(Internal_Intersections, ", ", External_Intersections)
        )
      )
    ) %>% separate_rows(to, sep = ",\\s*") %>% 
    mutate(
      type = ifelse(Internal_Intersections == to, "Internal", "External")
    ) %>% mutate(
      # below is to deal with the none cases in the network
      to = ifelse(is.na(type), from, to),
      type = ifelse(is.na(type), "None", type)
    ) %>% select(from, to, type)
  return(data_igraph)
}

igraphAssociations <- function(data){
  Associations(data) %>% graph_from_data_frame(data, directed = TRUE)
}

# Associations(read_xls("data.xls", sheet = 1))
# igraphAssociations(read_xls("data.xls", sheet = 1))

################################################################################
# Function to find associations of our associations to a given depth
# realistically, the max depth here is 2 for most cases

subAssociations <- function(data_igraph, rootID, depth) {
  depth_i <- 0
  horizon <- list(rootID) # what IDs constitute the root(s); impacted by depth
  observed <- c() # has this ID been seen already
  while (length(horizon) > 0 && depth_i < depth) {
    event_horizon <- c() # recursive horizon and trying to keep in theme
    for (i in horizon) {
      if (!i %in% observed) {
        observed <- observed %>% append(i)
        near <- neighbors(data_igraph, i, mode = "out") %>% as_ids() %>% setdiff(., i)
        event_horizon <- event_horizon %>% append(setdiff(near, observed))
      }
    }
    horizon <- unique(event_horizon)
    depth_i <- depth_i + 1
  }
  induced_subgraph(data_igraph, vids = observed)
}

# subAssociations(igraphAssociations(read_xls("data.xls", sheet = 1)), "DoARH-13Rxg", 3)
# Return node names: V(subAssociations(igraphAssociations(read_xls("data.xls", sheet = 1)), "DoARH-13Rxg", 3))$name
# Return TNAME (if there were ones in the data): V(subAssociations(igraphAssociations(read_xls("data.xls", sheet = 1)), "DoARH-13Rxg", 3))$TNAME
# Return internal/external intersections (edges, but without direction):
# - V(subAssociations(igraphAssociations(read_xls("data.xls", sheet = 1)), "DoARH-13Rxg", 3))$Internal_Intersections
# - V(subAssociations(igraphAssociations(read_xls("data.xls", sheet = 1)), "DoARH-13Rxg", 3))$External_Intersections
# Return geometry (multipoint sf coordinates): V(subAssociations(igraphAssociations(read_xls("data.xls", sheet = 1)), "DoARH-13Rxg", 3))$geometry
# Return type (use E bc it references the edges; ignore that that sounds messy and dumb and like it could've been done better so it doesn't reference edges):
# E(subAssociations(igraphAssociations(read_xls("data.xls", sheet = 1)), "DoARH-13Rxg", 3))$type

################################################################################
# Function to return property information for a root ID and it's associates.
# In the visnet, this will be truncated to only include some information
# returned by the function. But the other information returned can be useful
# in a table to go with the visnet.

getPropData <- function(data, dataPath, rootID, depth) {
  sub_nodes <- V(subAssociations(igraphAssociations(data), rootID, depth))$name
  for (i in 1:length(sub_nodes)) {
    ref_name <- ifelse(
      str_split_i(sub_nodes[[i]], "-", 1) %in% c("DoARH", "Grand List", "HDS", "IA", "JFO"),
      str_split_i(sub_nodes[[i]], "-", 1),
      ifelse(
        str_split_i(sub_nodes[[i]], "-", 1) == "GL", "Grand List",
        ifelse(
          str_split_i(sub_nodes[[i]], "-", 1) %in% c("DHCD", "VHCB", "VHFA"), "IA",
          "Columns not found."
        )
      )
    )
    # dataPath fixes problems in Shiny
    ref_data <- data.frame(ID = NA, addr = NA, SPAN = NA) %>% 
      plyr::rbind.fill(
        read_xls(dataPath, sheet = ref_name) %>% 
          filter(ID == sub_nodes[[i]]) %>% select(ID) %>% 
          cbind(
            data.frame(ID = NA, addr = NA, SPAN = NA) %>% 
              plyr::rbind.fill(
                read_xls(dataPath, sheet = ref_name) 
              ) %>% 
              filter(ID == sub_nodes[[i]]) %>% 
              mutate_if(is.character, str_to_title) %>% 
              mutate(
                SPAN = str_remove_all(SPAN, "-"),
                addr = str_replace_all(addr, " Vt", " VT")
              ) %>% select(-ID)
          ) 
      ) %>% filter(!is.na(ID))
    if (i == 1) {
      prop_data <- ref_data
    } else {
      prop_data <- prop_data %>% plyr::rbind.fill(ref_data)
    }
  }
  prop_data
}

simplePropData <- function(data, dataPath, rootID, depth) {
  data.frame(ID = NA, addr = NA, SPAN = NA) %>% 
    plyr::rbind.fill(
      getPropData(data, dataPath, rootID, depth)
    ) %>% select(ID, addr, SPAN) %>% unique() %>% 
    select(where(~!all(is.na(.x)))) %>% filter(!is.na(ID))
}

# test <- simplePropData(read_xls("data.xls", sheet = 1), "data.xls", "JFO-nTfoJ", 3)

# test[which(test$ID == "DoARH-13Rxg"),]
# test[which(test$ID == "GL-2MBOo"),]
# test[which(test$ID == "HDS-dKCRu"),]

################################################################################

toVisnet <- function(data, dataPath, rootID, depth) {
  sub_graph <- subAssociations(igraphAssociations(data), rootID, depth)
  sub_data <- simplePropData(data, dataPath, rootID, depth)
  labels <- as.data.frame(V(sub_graph)$name) %>% rename(to = 1) %>% filter(to != rootID) %>% 
    left_join(
      Associations(data) %>% filter(from == rootID) %>% select(to, type) %>% unique()
    )
  node_list <- tibble(
    id = V(sub_graph)$name,
    label = ifelse(V(sub_graph)$name == rootID, V(sub_graph)$name, NA),
    title = V(sub_graph)$name,
    # Not needed anymore; used for ifelse statement to identify ID type in tooltip
    #nodetype = ifelse(
    #  "None" %in% E(sub_graph)$name, "None",
    #  ifelse(
    #    substr(labels$to[labels$to %in% V(sub_graph)$name], 1, 2) == substr(rootID, 1, 2),
    #    "Internal", "External"
    #  )
    #),
    # paste(test %>% filter(ID == "HDS-dKCRu") %>% select(addr) %>% unlist(), collapse = "<br>")
    size = ifelse(V(sub_graph)$name == rootID, 50, 25),
    # labels$to[labels$to %in% V(sub_graph)$name]
  ) %>% rowwise() %>% 
    # addr_info: paste(prop_info %>% select(ID, addr) %>% unique() %>% unlist(), collapse = "<br>")
    # SPAN_info: paste(prop_info %>% select(ID, SPAN) %>% unique() %>% unlist(), collapse = "<br>")
    mutate(
      shape = ifelse(
        title == rootID, "circle",
        ifelse(
          substr(labels$to[labels$to %in% title], 1, 2) == substr(rootID, 1, 2),
          "square", "triangle"
        )
      ),
      color = ifelse(
        # VHFA green color for root
        title == rootID, "#378134",
        ifelse(
          substr(labels$to[labels$to %in% title], 1, 2) == substr(rootID, 1, 2),
          "#378134",
          ifelse(
            substr(labels$to[labels$to %in% title], 1, 2) == "GL", "goldenrod",
            ifelse(
              substr(labels$to[labels$to %in% title], 1, 2) == "HD", "darkorange",
              ifelse(
                substr(labels$to[labels$to %in% title], 1, 2) == "JF", "orchid",
                ifelse(
                  substr(labels$to[labels$to %in% title], 1, 2) %in% c("DHCD", "VHCB", "VHFA"), "slateblue",
                  ifelse(
                    substr(labels$to[labels$to %in% title], 1, 2) == "Do", "tomato",
                    "black"
                  )
                )
              )
            )
          )
        )
      ),
      title = paste0(
        "<p>ID: ", title, "</p>",
        "<p>Property Information: </p>",
        "<p>> Address", ifelse(length((sub_data[which(sub_data$ID %in% title),])$addr %>% unlist()) > 1, "es:", ":"), paste0("<br>&emsp;&emsp; | ", (sub_data[which(sub_data$ID %in% title),])$addr %>% unlist(), collapse = ""), "</p>"
      )
    )
  
  sub_graph_vis <- as_data_frame(sub_graph, what = "edges")
  has_associations <- sub_graph_vis %>% filter(type != "None")
  no_associations <- sub_graph_vis %>% filter(type == "None")
  sub_graph_visnet <- bind_rows(has_associations, no_associations) %>% 
    select(from, to) %>% mutate(smooth = TRUE)
  list(nodes = node_list, edges = sub_graph_visnet)
}

# this ID, GL-ULDLS, has the most external intersections (by number of characters in the cell string)
# print(toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "GL-ULDLS", 3)$nodes, n = 59)
# a smaller case: "DoARH-13Rxg"
# toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-13Rxg", 3)$edges 
# toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-13Rxg", 3)$nodes
# toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-13Rxg", 3)$nodes %>% select(ID) %>% unique()

# Case with no matches

#visNetwork(toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-1uosK", 3)$nodes,
#           toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-1uosK", 3)$edges) %>% 
#  visEdges(smooth = FALSE) %>% 
#  visOptions(
#    highlightNearest = TRUE, nodesIdSelection = TRUE
#  ) %>% 
#  visPhysics(stabilization = TRUE) %>%
#  visInteraction(
#    navigationButtons = TRUE, tooltipDelay = 0, tooltipStay = 1000
#  )

# Case with some matches

#visNetwork(toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-13Rxg", 3)$nodes,
#           toVisnet(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-13Rxg", 3)$edges) %>% 
#  visEdges(smooth = FALSE) %>% 
#  visOptions(
#    highlightNearest = TRUE, nodesIdSelection = list(
#      enabled = TRUE,
#      style = "width: 25%; height: 30px; background: #f8f8f8;
#      color: #378134; border: none; outline: none;" 
#    )
#  ) %>% 
#  visPhysics(stabilization = TRUE) %>%
#  visInteraction(
#    navigationButtons = TRUE, tooltipDelay = 0, tooltipStay = Inf
#  )

# It works! Yay! Now, let's make it a slightly smaller call when we put it in Shiny.

visVHFA <- function(data, dataPath, rootID, depth) {
  visNetwork(toVisnet(data, dataPath, rootID, depth)$nodes,
             toVisnet(data, dataPath, rootID, depth)$edges)%>% 
    visEdges(smooth = FALSE) %>% 
    visOptions(
      highlightNearest = TRUE, nodesIdSelection = list(
        enabled = FALSE,
        style = "width: 25%; height: 30px; background: #f8f8f8;
      color: #378134; border: none; outline: none;" 
      )
    ) %>% 
    visPhysics(stabilization = TRUE) %>%
    visInteraction(
      navigationButtons = TRUE, tooltipDelay = 0, tooltipStay = Inf
    )
}

# visVHFA(read_xls("data.xls", sheet = 1), "data.xls", "DoARH-13Rxg", 3)