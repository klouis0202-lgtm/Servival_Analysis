# ============================================================
# 中文導讀：對應 Ch2「Likelihood Construction, Inference for Parametric
# Survival Distributions」。下列頁碼以講義投影片右下角的 1–36 頁為準。
# 主線：右設限資料（第 3–8 頁）→ 分段指數模型（第 23–26 頁）
# → 估計精確度與區間選擇（第 27–30 頁）→ Monte Carlo 模擬延伸。
# 符號對照：U = 觀測時間；delta = 事件指標；r[j] = 區間事件數；
# W[j] = 區間總觀測風險時間；lambda.hat[j] = 該區間估計危險率。
# 本次只新增說明註解，保留所有原始可執行程式碼。
# 注意：下方 cuts.few 的第二次賦值會覆寫第一次，實際為 9 區間；
# cuts.many 為 7 區間，因此原有 few/many 名稱與圖例不代表實際多寡。
# Piecewise Exponential Approximation
# True hazard: h(t) = t^2
# ============================================================

# ------------------------------------------------------------
# 1. Generate right-censored survival data
# ------------------------------------------------------------

# 【第 3 頁：建立 (U_i, delta_i) 的模擬樣本】
# 固定亂數種子以便重現結果；n 是受試者數。
# censor.rate 是設限時間 C 的指數分布「速率」，不是設限比例 30%。
# 實際設限比例由 T 與 C 的分布共同決定，下面用樣本比例檢查。
set.seed(2026)

n <- 400
censor.rate <- 0.3

# True cumulative hazard:
# H(t) = t^3 / 3
#
# Since H(T) ~ Exp(1),
# T = (3E)^(1/3), E ~ Exp(1)

# 【第 24 頁：h、H、S 的關係；此處用於模擬資料】
# 真實危險率 h(t)=t^2，故 H(t)=積分_0^t s^2 ds=t^3/3，
# S(t)=exp(-H(t))=exp(-t^3/3)，f(t)=h(t)S(t)。
# rexp(n) 預設 rate=1；令 H(T)=E，即得 T=(3E)^(1/3)。
# 這是真實危險率隨時間上升的分布，並非單一固定速率的指數分布。
# 後續以區間常數近似 t^2，呼應第 22–23 頁弱結構模型的動機。
E <- rexp(n)
T <- (3 * E)^(1/3)

# Independent censoring time
# 【第 3、6–8 頁：獨立的隨機設限】
# 另行產生 C，讓 T 與 C 獨立，符合此講義的非資訊性設限假設。
# 生存分布與設限分布參數分離時，完整 likelihood 可分解；
# 推論生存分布時可使用乘積 f(U_i)^delta_i S(U_i)^(1-delta_i)。
C <- rexp(n, rate = censor.rate)

# Observed survival data
# 【第 3–5 頁：實際可觀測資料】
# pmin 逐人取最小值：U_i=min(T_i,C_i)。
# T_i<=C_i 時 delta_i=1，觀察到事件，likelihood 貢獻為 f(U_i)；
# 否則 delta_i=0，只知道 T_i>U_i，貢獻為 S(U_i)=1-F(U_i)。
# dat 只保存分析時可見的 U、delta；不是把設限當成事件。
U <- pmin(T, C)
delta <- as.integer(T <= C)

dat <- data.frame(
  U = U,
  delta = delta
)

head(dat)

# Proportion censored
mean(delta == 0)


# ------------------------------------------------------------
# 2. Function for fitting a piecewise exponential model
# ------------------------------------------------------------

# 【第 23–26 頁：給定切點，估計每一段的常數危險率】
# cuts 必須為正、遞增且不重複的內部切點；grid 為非負的評估時間。
# 函數假設 U、delta 等長，delta 為 0/1；此原始版本未額外檢查輸入。
# 模型有 J=length(cuts)+1 個參數，而非只有 cuts 的個數。
fit.pwe <- function(U, delta, cuts, grid) {

  # Partition:
  # [0, cut1), [cut1, cut2), ..., [last cut, Inf)

  # 【第 23–24 頁：切割時間軸】
  # 講義記為 v_0=0,...,v_{k+1}=Inf；R 從 1 開始編號，
  # 因此 R 的 v[j] 對應講義 v_{j-1}，v[j+1] 對應講義 v_j。
  # 最後一段延伸到無限大；並不因繪圖只到 t=2 就在 2 結束。
  v <- c(0, cuts, Inf)

  J <- length(v) - 1

  # 配置每段統計量：r 為事件個數、W 為總風險時間。
  # 第 25 頁的 likelihood 最後只需要這兩組統計量。
  r <- numeric(J)
  W <- numeric(J)

  for (j in 1:J) {

    # Number of observed events in interval j
    # 【第 25 頁】r_j=sum_i delta_i I(U_i 屬於 I_j)。
    # 只有事件且落在 [v[j],v[j+1]) 的個案才計入；設限不算事件。
    # 若時間恰等於切點，歸入右側區間。
    r[j] <- sum(
      delta == 1 &
      U >= v[j] &
      U < v[j + 1]
    )

    # Total observed time at risk in interval j
    # 【第 25 頁】W_j 是所有個案在本區間累積的觀測時間。
    # 每人貢獻 max(min(U_i,右端點)-左端點,0)：
    # 在左端前已結束追蹤 → 0；在區間內結束 → U_i-左端點；
    # 追蹤跨過右端點 → 整段長度。設限個案也貢獻設限前的時間。
    # 這與講義「區間內部分時間 + 跨過整段的時間」公式等價。
    W[j] <- sum(
      pmax(
        pmin(U, v[j + 1]) - v[j],
        0
      )
    )
  }

  # Piecewise hazard estimates
  # 【第 25–26 頁：最大概似估計】
  # L = 乘積_j lambda_j^r_j exp(-lambda_j W_j)，
  # log L = 總和_j [r_j log(lambda_j)-lambda_j W_j]。
  # 對 lambda_j 微分：r_j/lambda_j-W_j=0，故估計為 r_j/W_j。
  # 這也是第 10–11 頁單一指數分布估計式在各區間的應用。
  # W_j=0 表示無人提供觀測時間，無法估計，回傳 NA；
  # W_j>0 且 r_j=0 時回傳邊界值 0，此時不能用 log 常態近似。
  lambda.hat <- ifelse(
    W > 0,
    r / W,
    NA_real_
  )

  # Determine the interval corresponding to each grid point
  # 【第 24 頁：將每段參數還原成階梯危險率 h(t)】
  # findInterval 找出 grid 的每個時間屬於哪一段；
  # 隨後用 lambda.hat[index] 指派該段固定高度。
  # 這裡產生的是 h(t)，不是累積危險率 H(t) 或存活機率 S(t)。
  index <- findInterval(
    grid,
    c(0, cuts)
  )

  # Estimated hazard function over grid
  h.hat <- lambda.hat[index]

  # 回傳區間事件數、風險時間、參數估計，以及繪圖網格上的危險率。
  list(
    r = r,
    W = W,
    lambda.hat = lambda.hat,
    h.hat = h.hat
  )
}


# ------------------------------------------------------------
# 3. Define grid and user-specified cut points
# ------------------------------------------------------------

# 【第 23、27、30 頁：預先指定區間；grid 僅供繪圖】
# grid 的 200 點用來評估曲線，不代表模型有 200 個參數。
# 切點固定於資料之外，符合講義先指定區間的推論設定。
grid <- seq(
  0,
  2,
  length.out = 200
)

# Fewer intervals:
# 4 intervals = 4 hazard parameters
# 第一個版本：3 個切點、4 個區間；但稍後會被同名賦值覆寫。
cuts.few <- c(
  0.75,
  1.25,
  1.75
)


# Fewer intervals:
# 4 intervals = 4 hazard parameters
# 【原碼提醒】此處覆寫上面的 cuts.few，真正使用 8 個切點、9 段。
# 上方原始英文註解「4 intervals」只適用於第一次設定。
# 所以實際 fit.few 有 9 段，fit.many 有 7 段，圖例名稱相反。
# 本段註解保留這個差異，沒有擅自刪改原本的切點。
# grid 最後一點 t=2 會落入 [2,Inf)，採用第 9 段估計。
cuts.few <- c(
 0.25,
0.5, 
0.75,
1,
  1.25,
1.5,
  1.75,
2
)



# More intervals:
# 7 intervals = 7 hazard parameters
cuts.many <- c(
  0.50,
  0.75,
  1.00,
  1.25,
  1.50,
  1.75
)


# ------------------------------------------------------------
# 4. Fit both piecewise exponential models
# ------------------------------------------------------------

# 【第 25–26 頁】同一份觀測資料，用兩套固定切點分別計算 r_j/W_j。
# 比較差異來自切點設計；函數並沒有數值最佳化，因為有封閉解。
fit.few <- fit.pwe(
  U = U,
  delta = delta,
  cuts = cuts.few,
  grid = grid
)

fit.many <- fit.pwe(
  U = U,
  delta = delta,
  cuts = cuts.many,
  grid = grid
)


# ------------------------------------------------------------
# 5. Examine estimated parameters
# ------------------------------------------------------------

# 【第 26–28 頁：檢查估計表與精確度】
# 下列輸出對應第 28 頁表格的 r_j、W_j、lambda_hat_j。
# 第 26 頁觀測資訊矩陣為 diag(r_j/lambda_j^2)。
# 取逆並代入估計量，可估 Var(lambda_hat_j) 約為 r_j/W_j^2。
# Number of events in each interval
fit.few$r
fit.many$r

# Total time at risk in each interval
fit.few$W
fit.many$W

# Estimated piecewise hazards
fit.few$lambda.hat
fit.many$lambda.hat

# Approximate variance of log(lambda_hat_j)
# 【第 13、27 頁：delta method】
# Var(log(lambda_hat_j)) 約為 1/r_j，是「對數估計值的變異數」。
# 它不是 lambda_hat_j 本身的變異數，也不是標準差。
# r_j>0 且事件數足夠時，約 95% CI 可寫為
# lambda_hat_j * exp(正負 1.96/sqrt(r_j))（第 15 頁推廣）。
# r_j=0 時以下會顯示 Inf；這表示該 log 常態推論不能使用。
1 / fit.few$r
1 / fit.many$r


# ------------------------------------------------------------
# 6. Plot:
#    true hazard versus piecewise exponential estimates
# ------------------------------------------------------------

# 【第 24、29 頁：看分段危險率如何近似真實曲線】
# 實線為真值 t^2，兩條階梯線為單次樣本估計。
# range(...,na.rm=TRUE) 讓此圖的 y 軸涵蓋真值與兩組可用估計。
# type="s" 以階梯線呈現；網格有限，視覺跳點可能略有偏移。
# 第 30 頁另建議用 log 危險率觀察精確度；本碼仍畫原始尺度。
# 原圖例的 Fewer/More 請依上方覆寫提醒解讀。
ylim <- range(
  c(
    grid^2,
    fit.few$h.hat,
    fit.many$h.hat
  ),
  na.rm = TRUE
)

plot(
  grid,
  grid^2,
  type = "l",
  lwd = 2,
  ylim = ylim,
  xlab = "t",
  ylab = "Hazard"
)

lines(
  grid,
  fit.few$h.hat,
  type = "s",
  lty = 2,
  lwd = 2
)

lines(
  grid,
  fit.many$h.hat,
  type = "s",
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True hazard",
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(1, 2, 3),
  lwd = 2
)


# ============================================================
# 7. Repeated simulation
# ============================================================

# 【模擬延伸：呼應第 29–30 頁估計變異與區間選擇】
# 以下不是講義 AZT/ddI 實例，而是重複產生 n=400 的模擬樣本。
# B=1000 是獨立重複次數；每次重新生成 T、C，再估計兩套模型。
# 重新設同一種子，讓整個模擬可重現；第一次樣本會重現前面的樣本。
# est.few、est.many：每列是 grid 時間點，每欄是一次模擬。
# 每次兩模型共用同一批資料，方便比較不同切點的效果。
set.seed(2026)

B <- 1000

est.few <- matrix(
  NA_real_,
  nrow = length(grid),
  ncol = B
)

est.many <- matrix(
  NA_real_,
  nrow = length(grid),
  ncol = B
)

for (b in 1:B) {

  # Generate survival times from h(t) = t^2
  E.b <- rexp(n)

  T.b <- (3 * E.b)^(1/3)

  # Generate independent censoring times
  C.b <- rexp(
    n,
    rate = censor.rate
  )

  # Observed data
  U.b <- pmin(
    T.b,
    C.b
  )

  delta.b <- as.integer(
    T.b <= C.b
  )

  # Fit model with fewer intervals
  fit.few.b <- fit.pwe(
    U = U.b,
    delta = delta.b,
    cuts = cuts.few,
    grid = grid
  )

  # Fit model with more intervals
  fit.many.b <- fit.pwe(
    U = U.b,
    delta = delta.b,
    cuts = cuts.many,
    grid = grid
  )

  # 將第 b 次的整條估計曲線存入第 b 欄，供跨次數彙整。
  est.few[, b] <- fit.few.b$h.hat

  est.many[, b] <- fit.many.b$h.hat
}


# ------------------------------------------------------------
# 8. Monte Carlo mean estimated hazard
# ------------------------------------------------------------

# 【Monte Carlo 平均：估計方法的平均表現】
# 對每個時間點，將 B 次估計取平均，近似 E[h_hat(t)]。
# 平均曲線與 t^2 的差異可反映分段近似與有限樣本造成的偏差。
# 它不是把 B 份資料合併後重新估計，也不是單次資料的信賴區間。
# na.rm=TRUE 會排除無法估計的 NA，各時間點有效重複數可能不同。
# 此圖的 y 軸由初始的 t^2 決定，若平均估計超過範圍可能被截掉。
mean.few <- rowMeans(
  est.few,
  na.rm = TRUE
)

mean.many <- rowMeans(
  est.many,
  na.rm = TRUE
)

plot(
  grid,
  grid^2,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = "Hazard"
)

lines(
  grid,
  mean.few,
  type = "s",
  lty = 2,
  lwd = 2
)

lines(
  grid,
  mean.many,
  type = "s",
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True hazard",
    "Mean: fewer intervals",
    "Mean: more intervals"
  ),
  lty = c(1, 2, 3),
  lwd = 2
)


# ------------------------------------------------------------
# 9. Empirical standard deviation of estimated hazard
# ------------------------------------------------------------

# 【Monte Carlo 標準差：呼應第 26、29–30 頁精確度】
# apply(...,1,sd) 沿每一列取樣本標準差，衡量 h_hat(t) 跨樣本波動。
# 這是原始危險率尺度的經驗 SD，與前面的 1/r_j（log 尺度變異數）不同。
# 區間越細通常每段事件越少，變異可能越大，但此處兩組切點
# 並非單純的巢狀細分，不能保證每個時間點都有相同的大小順序。
sd.few <- apply(
  est.few,
  1,
  sd,
  na.rm = TRUE
)

sd.many <- apply(
  est.many,
  1,
  sd,
  na.rm = TRUE
)

plot(
  grid,
  sd.few,
  type = "l",
  lty = 2,
  lwd = 2,
  ylim = range(
    c(sd.few, sd.many),
    na.rm = TRUE
  ),
  xlab = "t",
  ylab = "Empirical SD of estimated hazard"
)

lines(
  grid,
  sd.many,
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(2, 3),
  lwd = 2
)


# ------------------------------------------------------------
# 10. Optional: empirical bias
# ------------------------------------------------------------

# 【Monte Carlo 偏差：講義分段近似概念的延伸】
# Bias_hat(t)=平均估計危險率-t^2；正值為平均高估，負值為低估。
# 水平 0 線代表無偏；這是模擬評估，真實資料通常不知道真實危險率。
true.hazard <- grid^2

bias.few <- mean.few - true.hazard

bias.many <- mean.many - true.hazard

plot(
  grid,
  bias.few,
  type = "l",
  lty = 2,
  lwd = 2,
  ylim = range(
    c(bias.few, bias.many),
    na.rm = TRUE
  ),
  xlab = "t",
  ylab = "Empirical Bias"
)

lines(
  grid,
  bias.many,
  lty = 3,
  lwd = 2
)

abline(
  h = 0,
  lty = 1
)

legend(
  "topleft",
  legend = c(
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(2, 3),
  lwd = 2
)


# ------------------------------------------------------------
# 11. Optional: empirical MSE
# ------------------------------------------------------------

# 【Monte Carlo 均方誤差：同時衡量偏差與變異】
# MSE_hat(t)=(1/B) sum_b [h_hat_b(t)-h(t)]^2（有 NA 時用有效次數）。
# R 對矩陣減去 true.hazard 時按欄重複此向量，所以每次模擬
# 都減去相同的逐時間真值；rowMeans 再對每列平均平方誤差。
# 理論 MSE=Var(h_hat)+Bias^2；MSE 越小代表該點總體誤差越小。
# 若有效次數為 m，這裡 MSE = Bias_hat^2 + (m-1)/m * 經驗SD^2，
# 因為 sd() 使用 m-1 作分母。細分區間不一定降低整體 MSE。
mse.few <- rowMeans(
  (est.few - true.hazard)^2,
  na.rm = TRUE
)

mse.many <- rowMeans(
  (est.many - true.hazard)^2,
  na.rm = TRUE
)

plot(
  grid,
  mse.few,
  type = "l",
  lty = 2,
  lwd = 2,
  ylim = range(
    c(mse.few, mse.many),
    na.rm = TRUE
  ),
  xlab = "t",
  ylab = "Empirical MSE"
)

lines(
  grid,
  mse.many,
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(2, 3),
  lwd = 2
)