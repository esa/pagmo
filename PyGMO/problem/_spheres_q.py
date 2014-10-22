def visualize(self, world_states):
    # If Vpython is not detected this method throws an excpetion
    try:
        import visual as v
    except ImportError:
        print(
            "error: No Visual Python Detected. This is needed for the visualization method")
        raise
    import numpy as np

    # We convert the input data in a numpy array
    world_states = np.array(world_states)

    # Initializes the display and the initial sphere positions
    scene = v.display(
        title="Spheres simulations",
        width=800,
        height=450,
        x=100, y=0,
        visible=True,
        autocenter=True,
        autoscale=False,
        exit=True,
    )

    # Creates the spheres
    bodies = [
        v.sphere(
            radius=0.1,
            color=v.color.white,
        )
        for _ in range(3)
    ]

    # Initialize spheres positionss
    for sph_id, sph in enumerate(bodies):
        sph.pos = world_states[0, sph_id * 3 + 1: (sph_id + 1) * 3 + 1]

    # Creates the trails
    trails = [v.curve(color=sc.color) for sc in bodies]
    for i in range(3):
        trails[i].append(pos=bodies[i].pos)

    dr = [
        np.sqrt((world_states[0, 1] - world_states[0, 4]) ** 2 + (world_states[0, 2] -
                world_states[0, 5]) ** 2 + (world_states[0, 3] - world_states[0, 6]) ** 2),
        np.sqrt((world_states[0, 1] - world_states[0, 7]) ** 2 + (world_states[0, 2] -
                world_states[0, 8]) ** 2 + (world_states[0, 2] - world_states[0, 9]) ** 2),
        np.sqrt((world_states[0, 4] - world_states[0, 7]) ** 2 + (world_states[0, 5] -
                world_states[0, 8]) ** 2 + (world_states[0, 6] - world_states[0, 9]) ** 2)
    ]

    # Creates the labels
    time_lbl = v.label(
        pos=(-9, -5, 0), text = "t = " + str(world_states[0, 0]))
    r_lbl = v.label(pos=(-
                         9, -
                         3, 0), text = "r12: " +
                    str(dr[0]) +
                    "\n" +
                    "r13: " +
                    str(dr[1]) +
                    "\n" +
                    "r23: " +
                    str(dr[2]))

    scene.autocenter = False

    # pause before starting animations
    v.rate(.5)

    #
    for i in range(1, world_states.shape[0]):

        for sph_id, sph in enumerate(bodies):
            sph.pos = world_states[i, sph_id * 3 + 1: (sph_id + 1) * 3 + 1]
        for j in range(3):
            trails[j].append(pos=bodies[j].pos, retain=100)

        dr = [
            np.sqrt((world_states[i, 1] - world_states[i, 4]) ** 2 + (world_states[i, 2] -
                    world_states[i, 5]) ** 2 + (world_states[i, 3] - world_states[i, 6]) ** 2),
            np.sqrt((world_states[i, 1] - world_states[i, 7]) ** 2 + (world_states[i, 2] -
                    world_states[i, 8]) ** 2 + (world_states[i, 3] - world_states[i, 9]) ** 2),
            np.sqrt((world_states[i, 4] - world_states[i, 7]) ** 2 + (world_states[i, 5] -
                    world_states[i, 8]) ** 2 + (world_states[i, 6] - world_states[i, 9]) ** 2)
        ]

        time_lbl.text = "t = " + str(world_states[i, 0])
        r_lbl.text = "r12: " + \
            str(dr[0]) + "\n" + "r13: " + \
            str(dr[1]) + "\n" + "r23: " + str(dr[2])
        v.rate(50)

    scene.visible = False

    for b in bodies:
        b.visible = False

    for t in trails:
        t.visible = False

    time_lbl.visible = False
    r_lbl.visible = False
    # del time_lbl
    # del r_lbl
    # del scene
