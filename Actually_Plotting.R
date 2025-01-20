ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% A_ResData +
  aes(color = I(FIA)) +
  scale_color_manual(
    name = "",
    breaks = c("FIA", "Not FIA"),
    labels = c("FIA", "Not FIA"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIA),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupA_FIA") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% A_FIB +
  aes(color = I(FIB)) +
  scale_color_manual(
    name = "",
    breaks = c("FIB", "No FIB"),
    labels = c("FIB", "No FIB"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIB),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupA_FIB") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(A_tree, layout = "circular", branch.length = 'none') %<+% A_FIC +
  aes(color = I(FIC)) +
  scale_color_manual(
    name = "",
    breaks = c("FIC", "No FIC"),
    labels = c("FIC", "No FIC"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIC),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupA_FIC") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(B1_tree, layout = "circular", branch.length = 'none') %<+% B1_FIA +
  aes(color = I(FIA)) +
  scale_color_manual(
    name = "",
    breaks = c("FIA", "No FIA"),
    labels = c("FIA", "No FIA"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIA),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupB1_FIA") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(B1_tree, layout = "circular", branch.length = 'none') %<+% B1_FIB +
  aes(color = I(FIB)) +
  scale_color_manual(
    name = "",
    breaks = c("FIB", "No FIB"),
    labels = c("FIB", "No FIB"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIB),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupB1_FIB") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(B1_tree, layout = "circular", branch.length = 'none') %<+% B1_FIC +
  aes(color = I(FIC)) +
  scale_color_manual(
    name = "",
    breaks = c("FIC", "No FIC"),
    labels = c("FIC", "No FIC"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIC),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupB1_FIC") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_FIA +
  aes(color = I(FIA)) +
  scale_color_manual(
    name = "",
    breaks = c("FIA", "No FIA", ""),
    labels = c("FIA", "No FIA", ""),
    values = c("blue", "white", "pink"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIA),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupB2_FIA") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_FIB +
  aes(color = I(FIB)) +
  scale_color_manual(
    name = "",
    breaks = c("FIB", "No FIB", ""),
    labels = c("FIB", "No FIB", ""),
    values = c("blue", "white", "pink"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIB),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupB2_FIB") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(B2_tree, layout = "circular", branch.length = 'none') %<+% B2_FIC +
  aes(color = I(FIC)) +
  scale_color_manual(
    name = "",
    breaks = c("FIC", "No FIC", ""),
    labels = c("FIC", "No FIC", "pink"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIC),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupB2_FIC") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(D_tree, layout = "circular", branch.length = 'none') %<+% D_FIA +
  aes(color = I(FIA)) +
  scale_color_manual(
    name = "",
    breaks = c("FIA", "Not FIA"),
    labels = c("FIA", "Not FIA"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIA),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupD_FIA") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(D_tree, layout = "circular", branch.length = 'none') %<+% D_FIB +
  aes(color = I(FIB)) +
  scale_color_manual(
    name = "",
    breaks = c("FIB", "No FIB"),
    labels = c("FIB", "No FIB"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIB),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupD_FIB") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )

ggtree(D_tree, layout = "circular", branch.length = 'none') %<+% D_FIC +
  aes(color = I(FIC)) +
  scale_color_manual(
    name = "",
    breaks = c("FIC", "No FIC"),
    labels = c("FIC", "No FIC"),
    values = c("blue", "white"),
    na.value = "grey"
  ) +
  geom_tippoint(
    mapping = aes(color = FIC),
    size = 1.5
  ) +
  geom_tiplab(
    color = 'black',
    offset = 1,
    size = 1,
    geom = "text",
    align = TRUE
  ) +
  ggtitle("PhylogroupD_FIC") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text=element_text(
      face = "bold",
      size = 10
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin()
  )
