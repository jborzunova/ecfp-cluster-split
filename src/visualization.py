from PIL import Image, ImageDraw
from PIL import ImageFont
from rdkit import Chem
from rdkit.Chem import Draw


def create_cluster_image(clusters, cluster_idx, compounds):
    try:
        grid_img = Draw.MolsToGridImage(
                                [Chem.MolFromSmiles(compounds[i][2]) for i in clusters[cluster_idx][:30]],
                                 legends=[compounds[i][1] for i in clusters[cluster_idx][:30]],
                                 molsPerRow=4, returnPNG=False, maxMols=160
                                )
    except Exception as e:
        print(f"Error generating grid image for cluster {cluster_idx}: {e}")
        return None

    # Add label at the top
    label_height = 40
    labeled_img = Image.new("RGB", (grid_img.width, grid_img.height + label_height), "white")
    draw = ImageDraw.Draw(labeled_img)
    # Попробуем загрузить системный жирный шрифт
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/tlwg/TlwgTypo-Bold.ttf", 18)
    except IOError:
        # Фолбэк, если шрифт не найден
        font = ImageFont.load_default()
        print("Bold font not found, using default font instead.")
    draw.text((10, 10), f"Cluster {cluster_idx}", fill="black", font=font)
    labeled_img.paste(grid_img, (0, label_height))

    return labeled_img


def stack_cluster_images(clusters, compounds):
    cluster_imgs = []

    for idx, _ in enumerate(clusters):
        img = create_cluster_image(clusters, idx, compounds)
        if img is not None:
            cluster_imgs.append(img)
        else:
            print(f"Cluster {idx} image generation failed or was empty.")

    if not cluster_imgs:
        raise ValueError("No valid cluster images to combine.")

    total_height = sum(img.height for img in cluster_imgs)
    max_width = max(img.width for img in cluster_imgs)

    combined = Image.new("RGB", (max_width, total_height), "white")

    y_offset = 0
    for img in cluster_imgs:
        combined.paste(img, (0, y_offset))
        y_offset += img.height
    return combined


def gen_and_save_single_page_pdf(clusters, compounds, filename="./pics/clustering.pdf"):
    img = stack_cluster_images(clusters, compounds)
    img.save(filename, "PDF", resolution=100.0)
    print(f'I saved {filename}')
