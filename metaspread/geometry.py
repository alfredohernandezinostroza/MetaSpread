"""Device geometry for organ-on-chip simulations (Phase 2).

A "wall mask" marks grid cells that agents (cancer and immune) cannot enter,
emulating the channels and chambers of a microfluidic device. The mask is a
boolean array the same spatial shape as the diffusible fields (2D or 3D); True
means an impassable wall.

Only consulted when ``enable_device_geometry`` is set, so default simulations
never build or read a mask and stay byte-identical.
"""
import numpy as np
import pandas as pd


def build_wall_mask(config, shape):
    """Return a boolean wall mask (True = impassable) for a field of ``shape``.

    Two sources, in priority order:

    - ``device_mask_path``: a file of 0/1 values (nonzero = wall) that must match
      ``shape`` exactly. ``.npy`` is loaded directly (any dimensionality); any
      other extension is read as a plain CSV grid with no header or index. This
      overrides the parametric channel.
    - otherwise a parametric straight **channel**: an open corridor runs along
      ``channel_axis``; ``channel_margin`` cells of wall line each side of the
      corridor on every other axis. ``channel_margin == 0`` leaves the domain
      fully open (a no-op mask), which is the inert default.
    """
    shape = tuple(int(s) for s in shape)

    path = getattr(config, "device_mask_path", "") or ""
    if path:
        if path.endswith(".npy"):
            mask = np.load(path)
        else:
            mask = pd.read_csv(path, header=None).to_numpy()
        mask = mask.astype(bool)
        if mask.shape != shape:
            raise ValueError(
                f"device_mask_path shape {mask.shape} does not match the field "
                f"shape {shape}."
            )
        return mask

    axis = int(getattr(config, "channel_axis", 0))
    margin = int(getattr(config, "channel_margin", 0))
    mask = np.zeros(shape, dtype=bool)
    if margin > 0:
        for other in range(len(shape)):
            if other == axis:
                continue  # the channel stays open along its own axis
            low = [slice(None)] * len(shape)
            low[other] = slice(0, margin)
            mask[tuple(low)] = True
            high = [slice(None)] * len(shape)
            high[other] = slice(shape[other] - margin, None)
            mask[tuple(high)] = True
    return mask
