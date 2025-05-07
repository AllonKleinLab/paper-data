import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Define function for cubic Bézier curve
def bezier_curve(t, P0, P1, P2, P3):
    output = ((1 - t) ** 3 * P0
              + 3 * (1 - t) ** 2 * t * P1
              + 3 * (1 - t) * t ** 2 * P2
              + t ** 3 * P3)
    return output


def ribbon_curves(x_start, x_end, y_start, y_end, w_start, w_end,
                  control_points_offset=3):
    # Upper Bezier curve
    P0 = [x_start, y_start]       # Start point (left)
    P1 = [x_start, y_start-control_points_offset] # Control point near the left
    P2 = [x_end, y_end+control_points_offset]   # Control point near the right
    P3 = [x_end, y_end]   # End point (right)
    t = np.linspace(0, 1, 400)
    bezier_x1 = bezier_curve(t, P0[0], P1[0], P2[0], P3[0])
    bezier_y = bezier_curve(t, P0[1], P1[1], P2[1], P3[1])
    # bezier1 = np.stack((bezier_x, bezier_y), axis=1)

    # Lower Bezier curve
    P0[0] = P0[0] + w_start       # Start point (left)
    P1[0] = P1[0] + w_start # Control point near the left
    P2[0] = P2[0] + w_end   # Control point near the right
    P3[0] = P3[0] + w_end   # End point (right)
    t = np.linspace(0, 1, 400)
    bezier_x2 = bezier_curve(t, P0[0], P1[0], P2[0], P3[0])

    return (bezier_y, bezier_x1, bezier_x2)


def plot_modified_sankey(top_rect_widths: dict, bot_rect_widths: dict,
                         ribbon_widths: list,
                         top_rect_unmatched=None, bot_rect_unmatched=None,
                         rect_height=0.5, scale_factor=50,
                         vertical_rect_spacing=2, horizontal_rect_spacing=None,
                         text_spacing=0.1, scalebar_size=50, figsize=(2, 1.5),
                         color_top='#8918d3', color_bot='#008c67',
                         color_ribbon='#888888'):
    """
    This function plots the modified sankey plot, used to show the number of
    genes shared between input cell types. Only one cell type is allowed in
    the top row.

    Parameters
    ----------
    top_rect_widths : dict, length must be 1
        blablabla description

    Returns
    -------
    f : figure object with sankey plot
    """

    # Create [top/bot]_rects_dict from [top/bot]_rect_widths
    top_rects_dict = {x:{'width':top_rect_widths[x]} for x in top_rect_widths}
    bot_rects_dict = {x:{'width':bot_rect_widths[x]} for x in bot_rect_widths}

    # Calculate additional plotting parameters
    rect_y_top = vertical_rect_spacing / 2
    rect_y_bot = -vertical_rect_spacing / 2
    control_points_offset = vertical_rect_spacing / 4

    if top_rect_unmatched is not None:
        top_widths = [
            (top_rects_dict[x]['width']/scale_factor
             + top_rect_unmatched[x]/scale_factor) for x in top_rects_dict]
    else:
        top_widths = [top_rects_dict[x]['width']/scale_factor for x in top_rects_dict]
    if bot_rect_unmatched is not None:
        bot_widths = [
            (bot_rects_dict[x]['width']/scale_factor
             + bot_rect_unmatched[x]/scale_factor) for x in bot_rects_dict]
    else:
        bot_widths = [bot_rects_dict[x]['width']/scale_factor for x in bot_rects_dict]
    if horizontal_rect_spacing is None:
        max_width = max(np.mean(top_widths), np.mean(bot_widths))
        horizontal_rect_spacing = 0.1 * max_width

    # ------------------------------------
    # Draw rectangles

    # Top rects
    leftmost_x_value_top = -0.5 * (
        np.sum(top_widths)
        + (len(top_widths) - 1) * horizontal_rect_spacing
    )
    current_x_top = leftmost_x_value_top
    for rect_label in top_rects_dict:
        # Plot rectangle
        rect_width = top_rects_dict[rect_label]['width'] / scale_factor
        top_rects_dict[rect_label]['starting_x'] = current_x_top
        top_rects_dict[rect_label]['current_x'] = current_x_top
        top_rects_dict[rect_label]['patch'] = patches.Rectangle(
            (current_x_top, rect_y_top),
            rect_width,
            rect_height,
            color=color_top,
            linewidth=0
        )
        if top_rect_unmatched is not None:
            unmatched_width = (top_rect_unmatched[rect_label]/scale_factor)
            top_rects_dict[rect_label]['patch2'] = patches.Rectangle(
                (current_x_top, rect_y_top),
                rect_width + unmatched_width,
                rect_height,
                color=color_top,
                linewidth=0,
                alpha=0.5,
            )
            rect_width += unmatched_width
        
        # Label rectangle
        text_y = rect_y_top + rect_height + text_spacing
        text_x = current_x_top + rect_width/2
        top_rects_dict[rect_label]['text_xy'] = (text_x, text_y)

        # Update current_x for next rectangle
        current_x_top += rect_width + horizontal_rect_spacing

    # Bottom rects
    leftmost_x_value_bot = -0.5 * (
        np.sum(bot_widths)
        + (len(bot_widths) - 1) * horizontal_rect_spacing
    )
    current_x_bot = leftmost_x_value_bot
    for rect_label in bot_rects_dict:
        rect_width = bot_rects_dict[rect_label]['width'] / scale_factor
        bot_rects_dict[rect_label]['starting_x'] = current_x_bot
        bot_rects_dict[rect_label]['current_x'] = current_x_bot
        bot_rects_dict[rect_label]['patch'] = patches.Rectangle(
            (current_x_bot, rect_y_bot),
            rect_width,
            -rect_height,
            color=color_bot,
            linewidth=0
        )
        if bot_rect_unmatched is not None:
            unmatched_width = (bot_rect_unmatched[rect_label]/scale_factor)
            bot_rects_dict[rect_label]['patch2'] = patches.Rectangle(
                (current_x_bot, rect_y_bot),
                rect_width + unmatched_width,
                -rect_height,
                color=color_bot,
                linewidth=0,
                alpha=0.5,
            )
            rect_width += unmatched_width
        
        # Label rectangle
        text_y = rect_y_bot - rect_height - text_spacing
        text_x = current_x_bot + rect_width/2
        bot_rects_dict[rect_label]['text_xy'] = (text_x, text_y)

        # Update current_x for next rectangle
        current_x_bot += rect_width + horizontal_rect_spacing

    # ------------------------------------
    # Draw ribbons

    ribbon_patches = []
    for ribbon in ribbon_widths:
        rect_top, rect_bot, w_top, w_bot, overlap_with_next = ribbon
        w_top /= scale_factor
        w_bot /= scale_factor
        overlap_with_next /= scale_factor
        current_x_top = top_rects_dict[rect_top]['current_x']
        current_x_bot = bot_rects_dict[rect_bot]['current_x']
        ribbon_patches.append(
            ribbon_curves(
                current_x_top, current_x_bot,
                rect_y_top, rect_y_bot,
                w_top, w_bot,
                control_points_offset=control_points_offset)
        )
        top_rects_dict[rect_top]['current_x'] += w_top - overlap_with_next
        bot_rects_dict[rect_bot]['current_x'] += w_bot

    # ------------------------------------
    # Draw scale bar

    scalebar_width = scalebar_size / scale_factor
    scalebar_height = -0.2 * rect_height
    scalebar_x = -scalebar_width / 2
    scalebar_y = rect_y_bot - 3.5*rect_height
    scalebar_patch = patches.Rectangle(
        (scalebar_x, scalebar_y),
        scalebar_width,
        scalebar_height,
        color='k',
        linewidth=0
    )
    
    # Label scalebar
    text_y = scalebar_y + scalebar_height - text_spacing
    text_x = scalebar_x + scalebar_width/2
    scalebar_text_xy = (text_x, text_y)

    # ------------------------------------
    # Plot

    f = plt.figure(figsize=figsize)
    ax = plt.subplot(1, 1, 1)

    for rect_label in top_rects_dict:
        ax.add_patch(top_rects_dict[rect_label]['patch'])
        if 'patch2' in top_rects_dict[rect_label]:
            ax.add_patch(top_rects_dict[rect_label]['patch2'])
        plt.text(top_rects_dict[rect_label]['text_xy'][0],
                top_rects_dict[rect_label]['text_xy'][1],
                rect_label,
                fontsize=8,
                color=color_top,
                horizontalalignment='center',
                verticalalignment='bottom')

    for rect_label in bot_rects_dict:
        ax.add_patch(bot_rects_dict[rect_label]['patch'])
        if 'patch2' in bot_rects_dict[rect_label]:
            ax.add_patch(bot_rects_dict[rect_label]['patch2'])
        plt.text(bot_rects_dict[rect_label]['text_xy'][0],
                bot_rects_dict[rect_label]['text_xy'][1],
                rect_label,
                color=color_bot,
                fontsize=8,
                horizontalalignment='center',
                verticalalignment='top')

    for ribbon in ribbon_patches:
        plt.fill_betweenx(ribbon[0], ribbon[1], ribbon[2],
                          color=color_ribbon, linewidth=0, alpha=0.5)
        
    ax.add_patch(scalebar_patch)
    plt.text(scalebar_text_xy[0],
             scalebar_text_xy[1],
             f'{scalebar_size} genes',
             fontsize=7,
             horizontalalignment='center',
             verticalalignment='top')

    plt.axis('off')
    plt.tight_layout()
    return f
