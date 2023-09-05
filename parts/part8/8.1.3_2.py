# 导入相关包
import torch
import torch.nn as nn
# 搭建LSTM模型
class LSTM(nn.Module):
    def __init__(self):
        super(LSTM, self).__init__()
        self.lstm = nn.LSTM(
            input_size=1,
            hidden_size=256,
            num_layers=2,
            batch_first=True
        )
        self.linear = nn.Linear(256, 1)
        self.activation= nn.ReLU()
    def forward(self, x):
        x, (ht, ct) = self.lstm(x)
        x = self.linear(x)
        x = self.activation(x)
        x = x.view(x.size(0), -1)
        return x
