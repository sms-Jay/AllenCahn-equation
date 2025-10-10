# -*- coding: gbk -*-
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

def load_data_from_file(filename):
    """���ļ���������"""
    with open(filename, 'r') as f:
        # ��ȡԪ����
        metadata = f.readline().split()
        time_steps = int(metadata[0])
        x_size = int(metadata[1])
        y_size = int(metadata[2])
        
        # ��ȡx����
        x_data = np.array([float(x) for x in f.readline().split()])/50-1.0
        
        # ��ȡy����
        y_data = np.array([float(y) for y in f.readline().split()])/50-1.0
        
        # ��ʼ��u����
        u_data = np.zeros((time_steps, x_size, y_size))
        
        # ��ȡu����
        current_t = -1
        current_x = 0
        
        for line in f:
            line = line.strip()
            if line.startswith('t='):
                # ��ʱ�䲽
                current_t = int(line[2:])
                current_x = 0
            elif line:
                # ������
                values = [float(v) for v in line.split()]
                u_data[current_t, current_x, :] = values
                current_x += 1
    
    return u_data, x_data, y_data

def create_heatmaps(u, x, y, time_indices=None, time=0):
    """������ͼ"""
    if time_indices is None:
        time_indices = range(min(5, u.shape[0]))  # Ĭ����ʾǰ5��ʱ�䲽
    
    n_plots = len(time_indices)
    
    # ������ͼ
    fig, axes = plt.subplots(1, n_plots, figsize=(5*n_plots, 4))
    
    if n_plots == 1:
        axes = [axes]
    
    # ȷ����ɫ��Χ������ͼʹ����ͬ����ɫ��Χ��
    vmin = np.min(u)
    vmax = np.max(u)
    
    # ����ÿ��ʱ�䲽����ͼ
    images = []
    for i, t in enumerate(time_indices):
        im = axes[i].imshow(u[t], extent=[x[0], x[-1], y[0], y[-1]], 
                           origin='lower', cmap='viridis', aspect='auto',
                           vmin=vmin, vmax=vmax)
        axes[i].set_title(f't = {t}')
        axes[i].set_xlabel('x')
        axes[i].set_ylabel('y')
        images.append(im)
    
    # �����ɫ��
    plt.colorbar(images[0], ax=axes, shrink=0.6)
    
    plt.tight_layout()
    plt.savefig(f"heatmaps_t_{time}.png", dpi=1000)



def create_contour_plots(u, x, y, time_indices=None,time=1):
    """�����ȸ���ͼ"""
    if time_indices is None:
        time_indices = range(min(5, u.shape[0]))  # Ĭ����ʾǰ5��ʱ�䲽
    
    n_plots = len(time_indices)
    
    # ������ͼ
    fig, axes = plt.subplots(1, n_plots, figsize=(5*n_plots, 4))
    
    if n_plots == 1:
        axes = [axes]
    
    # ��������
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # ȷ����ɫ��Χ������ͼʹ����ͬ����ɫ��Χ��
    vmin = np.min(u)
    vmax = np.max(u)
    
    # ����ÿ��ʱ�䲽�ĵȸ���ͼ
    for i, t in enumerate(time_indices):
        contour = axes[i].contourf(X, Y, u[t], 20, cmap='viridis', vmin=vmin, vmax=vmax)
        axes[i].set_title(f't = {t}')
        axes[i].set_xlabel('x')
        axes[i].set_ylabel('y')
        
        # �����ɫ��
        plt.colorbar(contour, ax=axes[i])
    
    plt.tight_layout()
    plt.savefig(f"contour_plots_t_{time}.png", dpi=1000)

def main():
    """������"""
    # ��������
    data_file = "H:/undergraduate/scientific research/Allen Cahn equation/Allen_CahnSolver/u_test_data.txt"
    
    if not os.path.exists(data_file):
        print(f"�����ļ� {data_file} ������")
        print("��������C++������������")
        return
    
    u, x, y = load_data_from_file(data_file)
    print(f"������״: u({u.shape}), x({x.shape}), y({y.shape})")
    
    # ѡ��Ҫ���ӻ���ʱ���
    for time in [0,10,20,30,40,50,60,70,80,90,100]:

        time_indices = [time]
        create_heatmaps(u, x, y, time_indices,time)
    
        create_contour_plots(u, x, y, time_indices,time)

if __name__ == "__main__":
    main()