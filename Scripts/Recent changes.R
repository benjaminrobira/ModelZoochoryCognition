What has changed since last sving (CEFE version):

- I have now added a new index (alignment) to measure routes!
- I have now added that dispersal only occur if space is "free". Free space is defined as no tree in a circular buffer whose radius is defined so that the area of all circular buffers around the trees equal that of the map.
- I have removed learning, and new seedling are "learned" only if tree that dies is a known one.
- I have imposed that the seedling could only produce the next year
- I have now added movement from tree to encountered tree while planning to go to the target (instead of linear to the target with telescopic arms)
