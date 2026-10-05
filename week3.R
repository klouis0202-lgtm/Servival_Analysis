#假設T~Exp(lambda) ; C~Exp(mu)
# 上面是說明:存活時間 T 服從參數為 lambda 的指數分布,設限時間 C 服從參數為 mu 的指數分布

lambda = 0.5        # 設定真實的存活時間指數分布參數 lambda(事件發生率)
mu = 0.3             # 設定真實的設限時間指數分布參數 mu(設限發生率)
#------------------------------------------------------------------------------
B = 2000             # 模擬重複次數(replication),要模擬 2000 次
n = 200               # 每次模擬的樣本數(受試者人數)
Lambda = rep(NA, B)         # 建立長度為 B 的向量,存放每次模擬得到的 MLE 估計值 λ̂ = r/W,初始值皆為 NA
Lambda_Naive = rep(NA, B)   # 建立長度為 B 的向量,存放每次模擬得到的「naive」估計值 n/W,初始值皆為 NA

for(i in 1:B){                       # 開始迴圈,重複模擬 B 次(i 從 1 到 2000)
  
  T = rexp(n, rate = lambda)          # 產生 n 個服從 Exp(lambda) 的真實存活時間 T(未觀測)
  C = rexp(n, rate = mu)              # 產生 n 個服從 Exp(mu) 的潛在設限時間 C(未觀測)
  # observed data
  U = pmin(T, C)                      # 觀測時間 U = min(T, C),取 T 與 C 中較小者(對應投影片的 Uᵢ = min{Tᵢ, Cᵢ})
  delta = as.integer(T<=C)            # 設限指示變數 δ:若 T ≤ C(事件先發生)則為 1,否則(被設限)為 0
  #cbind(C, T, U, delta)              # (被註解掉)原本可用來檢視 C, T, U, delta 四欄並排比較
  
  dat = data.frame(U, delta = delta)  # 把觀測到的 U 與 δ 組成資料框,模擬真實情況下我們只能看到的資料
  
  #mean(delta == 0) # censoring rate  # (被註解掉)原本用來計算此次模擬的設限比例(δ=0 的比例)
  
  r = sum(delta)                      # r = 未設限(事件)個數,即 Σδᵢ,對應投影片的 r
  W = sum(U)                          # W = 總觀察時間,即 ΣUᵢ,對應投影片的 W
  Lambda[i] = lambda.hat = r/W        # 計算 MLE:λ̂ = r/W,同時存入 lambda.hat 這個變數,並存進 Lambda 向量的第 i 個位置
  Lambda_Naive[i] = n/W               # 計算「naive」估計值:用總人數 n(而非事件數 r)除以 W,存進 Lambda_Naive 向量
  #Naive estimation
  lambda.naive = n/W                  # 重複計算一次 naive 估計值(僅存成一個暫存變數,每次迴圈會被覆蓋,不影響結果)
}                                     # 結束迴圈

######################################
hist(Lambda)          # 畫出 2000 次模擬所得 MLE 估計值 λ̂ 的直方圖,觀察其抽樣分布
hist(Lambda_Naive)    # 畫出 2000 次模擬所得 naive 估計值的直方圖,與 MLE 分布做比較

var(Lambda)           # 計算 2000 次模擬中 MLE 估計值 λ̂ 的變異數,評估估計量的變異程度
