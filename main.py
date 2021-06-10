# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy import stats

path_root = r".\Sheet"
bs_sheet = pd.read_csv(path_root + "/BS_Ant_Consolidate.csv", encoding="gbk", index_col=0, na_values='Nan').T
cf_sheet = pd.read_csv(path_root + "/CF_Ant_Consolidate.csv", encoding="gbk", index_col=0, na_values='Nan').T
is_sheet = pd.read_csv(path_root + "/IS_Ant_Consolidate.csv", encoding="gbk", index_col=0, na_values='Nan').T
notes_sheet = pd.read_csv(path_root + "/Notes_Ant_Consolidate.csv", encoding="gbk", index_col=0, na_values='Nan').T


class Quality:
    def __init__(self):
        self.DSR = 0.0
        # days sales receivable index =\
        # receivables(t)/sales(t) / receivables(t-1)/sales(t-1)
        self.GMI = 0.0
        # gross margin index =\
        # gross margin(t-1) / gross margin(t)
        self.AQI = 0.0
        # asset quality index =\
        # (1-(PPE(t)+CA(t))/TA(t)) / (1-(PPE(t-1)+CA(t-1))/TA(t-1))
        self.SGI = 0.0
        # sales growth index = \
        # sales(t) / sale(t-1)
        self.DEPI = 0.0
        # depreciation index =\
        # depreciation rate(t-1) / depreciation rate(t)
        # depreciation rate = depreciation/(depreciation+PPE)
        self.SGAI = 0.0
        # sale,general,and administrative expenses index =\
        # (SGA(t)/sales(t)) / (SGA(t-1)/sales(t-1))
        self.Accurals = 0.0
        # =\
        # (income before extraordinary items - cfo) / total assets
        self.LEVI = 0.0
        # leverage index =\
        # leverage(t) / leverage(t-1)
        # leverage ratio = D/A
        # 疑问 是D/A还是L/A D是付息负债 L是总负债

        self.A = 0.0
        # = WC / TA
        self.B = 0.0
        # = RE / TA
        self.C = 0.0
        # = EBIT / TA
        self.D = 0.0
        # = MV of Equity / BV of Debt
        self.E = 0.0
        # = Revenue / TA

        self.current = ""
        self.last = ""

        self.placeholder = "======"
    # 分类思路: 复合指标直接在类中作为公共变量, 而初级指标直接使用函数
    # 跨期指标使用set_period方法, 单期指标在函数参数中使用period

    def set_period(self):
        if self.current == "2020.6.30":
            self.last = "2019.6.30"
        elif self.current == "2019.12.31":
            self.last = "2018.12.31"
        elif self.current == "2018.12.31":
            self.last = "2017.12.31"
        else:
            self.current = "2020.6.30"
            self.last = "2019.12.31"

    def set_period_custom(self, current, last):
        self.current = current
        self.last = last

    def debt(self, period):
        debt = bs_sheet["短期借款"][period]+bs_sheet["长期借款"][period] +\
               bs_sheet["其他应付款"][period]+bs_sheet["租赁负债"][period]
        return debt

    def ebit(self, period):
        # 用 营业利润+财务费用 的方法拟合 EBIT
        ebit = is_sheet["营业利润"][period] + is_sheet["财务费用"][period]
        return ebit

    def tax_burden(self, period):
        # tax burden = NI / EBT
        tax_burden = is_sheet["净利润"][period]/is_sheet["利润总额"][period]
        return tax_burden

    def interest_burden(self, period):
        # interest burden = EBT / EBIT
        interest_burden = is_sheet["利润总额"][period]/self.ebit(period)
        return interest_burden

    def ebit_margin(self, period):
        # EBIT margin = EBIT / R
        ebit_margin = self.ebit(period)/is_sheet["营业收入"][period]
        return ebit_margin

    def asset_turnover(self, period):
        # asset turnover = R / A
        asset_turnover = is_sheet["营业收入"][period]/bs_sheet["资产总计"][period]
        return asset_turnover

    def financial_leverage(self, period):
        # financial leverage = A / E
        financial_leverage = bs_sheet["资产总计"][period]/bs_sheet["所有者权益合计"][period]
        return financial_leverage

    def roe(self, period):
        # ROE = NI / E
        roe = is_sheet["净利润"][period]/bs_sheet["所有者权益合计"][period]
        return roe

    def noa(self, period):
        # net operating assets(NOA) = (total assets - cash)-(total liabilities - total debt)
        noa = (bs_sheet["资产总计"][period]-bs_sheet["货币资金"][period]) - (bs_sheet["负债总计"][period]-self.debt(period))
        return noa

    def arbs(self):
        # aggregated accruals under b/s approach
        # accrual ratio
        if self.current != self.last:
            aabs = self.noa(self.current) - self.noa(self.last)
            arbs = aabs/((self.noa(self.current) + self.noa(self.last))/2)
        else:
            arbs = np.nan
        return arbs

    def arcf(self):
        # aggregated accruals under c/f approach
        # accrual ratio
        if self.current != self.last:
            aacf = is_sheet["净利润"][self.current] - (cf_sheet["经营活动产生/(使用)的现金流量净额"][self.current] +
                                                    cf_sheet["投资活动使用的现金流量净额"][self.current])
            arcf = aacf/((self.noa(self.current) + self.noa(self.last))/2)
        else:
            arcf = np.nan
        return arcf

    def m_score_ready(self):
        if self.current != self.last:
            self.DSR = (bs_sheet["应收账款"][self.current]/is_sheet["营业收入"][self.current]) /\
                (bs_sheet["应收账款"][self.last]/is_sheet["营业收入"][self.last])
            self.GMI = self.gross_profit(self.last)/is_sheet["营业收入"][self.last] /\
                       self.gross_profit(self.current)/is_sheet["营业收入"][self.current]
            self.AQI = (1.0-(bs_sheet["固定资产"][self.current]+bs_sheet["流动资产合计"][self.current])/bs_sheet["资产总计"][self.current]) /\
                       (1.0-(bs_sheet["固定资产"][self.last]+bs_sheet["流动资产合计"][self.last])/bs_sheet["资产总计"][self.last])
            self.SGI = is_sheet["营业收入"][self.current] / is_sheet["营业收入"][self.last]
            self.DEPI = notes_sheet["累计折旧"][self.last]/(notes_sheet["累计折旧"][self.last]+bs_sheet["固定资产"][self.last]) /\
                        notes_sheet["累计折旧"][self.current]/(notes_sheet["累计折旧"][self.current]+bs_sheet["固定资产"][self.current])
            self.SGAI = (is_sheet["销售费用"][self.current]+is_sheet["管理费用"][self.current])/is_sheet["营业收入"][self.current] /\
                        (is_sheet["销售费用"][self.last]+is_sheet["管理费用"][self.last])/is_sheet["营业收入"][self.last]
            self.Accurals = (is_sheet["利润总额"][self.current] - notes_sheet["非经常性损益净额"][self.current] -
                             cf_sheet["经营活动产生/(使用)的现金流量净额"][self.current]) / bs_sheet["资产总计"][self.current]
            self.LEVI = self.debt(self.current)/bs_sheet["资产总计"][self.current] / self.debt(self.last)/bs_sheet["资产总计"][self.last]
        else:
            self.DSR = np.nan
            self.GMI = np.nan
            self.AQI = np.nan
            self.SGI = np.nan
            self.DEPI = np.nan
            self.SGAI = np.nan
            self.Accurals = np.nan
            self.LEVI = np.nan
        # self.LEVI = (bs_sheet["短期借款"][self.current]+bs_sheet["长期借款"][self.current]+bs_sheet["其他应付款"][self.current]+bs_sheet["租赁负债"][self.current])/bs_sheet["资产总计"][self.current] /\
        #             (bs_sheet["短期借款"][self.last]+bs_sheet["长期借款"][self.last]+bs_sheet["其他应付款"][self.last]+bs_sheet["租赁负债"][self.last])/bs_sheet["资产总计"][self.last]

    def z_score_ready(self):
        self.A = self.working_capital(self.current)/bs_sheet["资产总计"][self.current]
        self.B = bs_sheet["盈余公积"][self.current]/bs_sheet["资产总计"][self.current]
        self.C = self.ebit(self.current)/bs_sheet["资产总计"][self.current]
        self.D = bs_sheet["所有者权益合计"][self.current]/bs_sheet["负债总计"][self.current]
        self.E = is_sheet["营业收入"][self.current]/bs_sheet["资产总计"][self.current]

    def m_score(self):
        # Beneish Model
        # 查单侧正态分布表 M分越高越糟糕
        self.m_score_ready()
        m_score = -4.84 + 0.92 * self.DSR + 0.528 * self.GMI + 0.404 * self.AQI + 0.892 * self.SGI\
                  + 0.115 * self.DEPI - 0.172 * self.SGAI + 4.67 * self.Accurals - 0.327 * self.LEVI
        return m_score

    def z_score(self):
        # Altman:Bankruptcy Prediction Models
        # 小于1.81为糟糕 大于3.00为优秀
        self.z_score_ready()
        z_score = 1.2 * self.A + 1.4 * self.B + 3.3 * self.C + 0.6 * self.D + 1.0 * self.E
        return z_score

    def quick_asset(self,period):
        quick_asset = bs_sheet["流动资产合计"][period] - bs_sheet["其他流动资产"][period] - bs_sheet["一年内到期的非流动资产"][period] -\
                      bs_sheet["买入返售金融资产"][period]
        return quick_asset

    def working_capital(self, period):
        working_capital = bs_sheet["流动资产合计"][period] - bs_sheet["流动负债合计"][period]
        return working_capital

    def inventory(self, period):
        inventory = 0.0
        return inventory

    def gross_profit(self, period):
        gross_profit = is_sheet["营业收入"][period]-is_sheet["减：营业成本"][period]
        return gross_profit

    def roa(self,period):
        roa = self.ebit(period) / bs_sheet["资产总计"][period]
        return roa

    def total_expense(self,period):
        total_expense = is_sheet["减：营业成本"][period]+is_sheet["销售费用"][period]+is_sheet["管理费用"][period]+\
                        is_sheet["研发费用"][period]+is_sheet["财务费用"][period]
        return total_expense

    def RPCE(self,period):
        rpce = is_sheet["营业利润"][period]/self.total_expense(period)
        return rpce

    def average_RA(self):
        if self.current != self.last:
            average_RA = (bs_sheet["应收账款"][self.current] + bs_sheet["应收账款"][self.last]) / 2
        else:
            average_RA = np.nan
        return average_RA

    def average_CA(self):
        if self.current != self.last:
            average_CA = (bs_sheet["流动资产合计"][self.current] + bs_sheet["流动资产合计"][self.last]) / 2
        else:
            average_CA = np.nan
        return average_CA

    def average_asset(self):
        if self.current != self.last:
            average_asset = (bs_sheet["资产总计"][self.current] + bs_sheet["资产总计"][self.last]) / 2
        else:
            average_asset = np.nan
        return average_asset

    def average_equity(self):
        if self.current != self.last:
            average_equity = (bs_sheet["所有者权益合计"][self.current] + notes_sheet["期初所有者权益"][self.current]) / 2
        else:
            average_equity = np.nan
        return average_equity

    def result_sheet(self):
        result_list = []
        period_list = ["2020.6.30", "2019.12.31", "2018.12.31", "2017.12.31"]
        period_add = ["2020.6.30", "2019.12.31", "2019.6.30"]
        for i in range(len(period_list)):
            self.current = period_list[i]
            if i == len(period_list)-1:
                self.last = period_list[i]
                index_name = ["2016.12.31"+'-'+self.current]
            else:
                self.last = period_list[i+1]
                index_name = [self.last+'-'+self.current]
            m_score = self.m_score()
            z_score = self.z_score()
            result_dic = {
                "Quality Analysis": self.placeholder,
                "Beneish Model": np.nan,
                "DSR": self.DSR,
                "GMI": self.GMI,
                "AQI": self.AQI,
                "SGI": self.SGI,
                "DEPI": self.DEPI,
                "SGAI": self.SGAI,
                "Accruals": self.Accurals,
                "LEVI": self.LEVI,
                "M_score": m_score,
                "M_Prob": stats.norm.cdf(m_score),
                "Altman Model": np.nan,
                "A:WC/TA": self.A,
                "B:RE/TA": self.B,
                "C:EBIT/TA": self.C,
                "D:E/D": self.D,
                "E:Revenue/TA": self.E,
                "Z_score": z_score,
                "ARBS": self.arbs(),
                "ARCF": self.arcf(),
                "Debt-Paying ability": self.placeholder,
                "Current Asset": bs_sheet["流动资产合计"][self.current],
                "Current Liability": bs_sheet["流动负债合计"][self.current],
                "Current Ratio": bs_sheet["流动资产合计"][self.current]/bs_sheet["流动负债合计"][self.current],
                "Inventory": self.inventory(self.current),
                "Quick Asset1": bs_sheet["流动资产合计"][self.current] - self.inventory(self.current),
                "Quick Ratio1": (bs_sheet["流动资产合计"][self.current] - self.inventory(self.current))/bs_sheet["流动负债合计"][self.current],
                "Quick Asset2": self.quick_asset(self.current),
                "Quick Ratio2": self.quick_asset(self.current)/bs_sheet["流动负债合计"][self.current],
                "Long-term Liability": bs_sheet["非流动负债合计"][self.current],
                "Working Capital": self.working_capital(self.current),
                "LL/WC": bs_sheet["非流动负债合计"][self.current]/self.working_capital(self.current),
                "EBIT": self.ebit(self.current),
                "I Protection Multiples": self.ebit(self.current)/is_sheet["其中：利息费用"][self.current],
                "Liability": bs_sheet["负债总计"][self.current],
                "Asset": bs_sheet["资产总计"][self.current],
                "Debt-to-Asset Ratio": bs_sheet["负债总计"][self.current]/bs_sheet["资产总计"][self.current],
                "Profitability Analysis": self.placeholder,
                "Operating Revenue": is_sheet["营业收入"][self.current],
                "Operating Cost": is_sheet["减：营业成本"][self.current],
                "Gross Profit": self.gross_profit(self.current),
                "Gross Profit Margin": self.gross_profit(self.current)/is_sheet["营业收入"][self.current],
                "Net Profit": is_sheet["净利润"][self.current],
                "Net Profit Margin": is_sheet["净利润"][self.current]/is_sheet["营业收入"][self.current],
                "ROA": self.roa(self.current),
                "RPCE": self.RPCE(self.current),
                "DuPont Analysis": np.nan,
                "ROE": self.roe(self.current),
                "Tax Burden": self.tax_burden(self.current),
                "Interest Burden": self.interest_burden(self.current),
                "EBIT Margin": self.ebit_margin(self.current),
                "Asset Turnover": self.asset_turnover(self.current),
                "Financial Leverage": self.financial_leverage(self.current),
                "Operating Analysis": self.placeholder,
                "Average Receivables": self.average_RA(),
                "RA Turnover Ratio": is_sheet["营业收入"][self.current]/self.average_RA(),
                "RA Turnover Days": 365.0 / is_sheet["营业收入"][self.current]/self.average_RA(),
                "Average Currents": self.average_CA(),
                "CA Turnover Ratio": is_sheet["营业收入"][self.current]/self.average_CA(),
                "CA Turnover Days": 365.0 / is_sheet["营业收入"][self.current]/self.average_CA()
                          }
            result_cache = pd.DataFrame(result_dic, index=index_name)
            result_list.append(result_cache)
        result_sheet = pd.concat(result_list)
        result_sheet.T.to_csv(path_or_buf="./result.csv")
        return result_sheet

    def bs_structure_sheet(self):
        bs_sheet_copy = bs_sheet.copy(deep=True)
        bs_structure_A = bs_sheet_copy.T.iloc[:33, ]
        # bs_structure_L = bs_sheet_copy.T.iloc[34:60, ]
        # bs_structure_E = bs_sheet_copy.T.iloc[61:71, ]
        bs_structure_LE = bs_sheet_copy.T.iloc[34:72, ]
        for i in bs_sheet_copy.index:
            bs_structure_A[i] = bs_structure_A[i]/bs_structure_A[i][-1]
            bs_structure_LE[i] = bs_structure_LE[i]/bs_structure_LE[i][-1]
            # bs_structure_L[period_list[i]] = bs_structure_L[period_list[i]]/bs_structure_L[period_list[i]][-1]
            # bs_structure_E[period_list[i]] = bs_structure_E[period_list[i]]/bs_structure_E[period_list[i]][-1]
        # bs_structure_sheet = pd.concat([bs_structure_A, bs_structure_L, bs_structure_E])
        bs_structure_sheet = pd.concat([bs_structure_A, bs_structure_LE])
        bs_structure_sheet.to_csv(path_or_buf="./bs_structure_sheet.csv",encoding="gbk")
        return bs_structure_sheet

    def is_structure_sheet(self):
        is_sheet_copy = is_sheet.copy(deep=True)
        is_structure_1 = is_sheet_copy.T.iloc[:24, ]
        for i in is_sheet_copy.index:
            is_structure_1[i] = is_structure_1[i]/is_structure_1[i][0]
        is_structure_sheet = pd.concat([is_structure_1])
        is_structure_sheet.to_csv(path_or_buf="./is_structure_sheet.csv",encoding="gbk")
        return is_structure_sheet

    def cf_structure_sheet(self):
        cf_sheet_copy = cf_sheet.copy(deep=True)
        cf_structure_cfo_in = cf_sheet_copy.T.iloc[:4, ]
        cf_structure_cfo_out = cf_sheet_copy.T.iloc[4:11, ]
        cf_structure_cfi_in = cf_sheet_copy.T.iloc[12:19, ]
        cf_structure_cfi_out = cf_sheet_copy.T.iloc[19:24, ]
        cf_structure_cff_in = cf_sheet_copy.T.iloc[25:31, ]
        cf_structure_cff_out = cf_sheet_copy.T.iloc[32:36, ]
        cf_structure_list = [cf_structure_cfo_in,cf_structure_cfo_out,cf_structure_cfi_in,
                             cf_structure_cfi_out,cf_structure_cff_in,cf_structure_cff_out]
        for i in cf_sheet.index:
            cf_structure_cfo_in[i] = cf_structure_cfo_in[i]/cf_structure_cfo_in[i][-1]
            cf_structure_cfo_out[i] = cf_structure_cfo_out[i]/cf_structure_cfo_out[i][-1]
            cf_structure_cfi_in[i] = cf_structure_cfi_in[i]/cf_structure_cfi_in[i][-1]
            cf_structure_cfi_out[i] = cf_structure_cfi_out[i]/cf_structure_cfi_out[i][-1]
            cf_structure_cff_in[i] = cf_structure_cff_in[i]/cf_structure_cff_in[i][-1]
            cf_structure_cff_out[i] = cf_structure_cff_out[i]/cf_structure_cff_out[i][-1]
        cf_structure_sheet = pd.concat(cf_structure_list)
        cf_structure_sheet.loc["总计"] = np.nan
        cf_structure_sheet.loc["现金流入总计"] = cf_sheet["经营活动现金流入小计"] + cf_sheet["投资活动现金流入小计"] + cf_sheet["筹资活动现金流入小计"]
        cf_structure_sheet.loc["CFO_IN"] = cf_sheet["经营活动现金流入小计"] / cf_structure_sheet.loc["现金流入总计"]
        cf_structure_sheet.loc["CFI_IN"] = cf_sheet["投资活动现金流入小计"] / cf_structure_sheet.loc["现金流入总计"]
        cf_structure_sheet.loc["CFF_IN"] = cf_sheet["筹资活动现金流入小计"] / cf_structure_sheet.loc["现金流入总计"]
        cf_structure_sheet.loc["现金流出总计"] = cf_sheet["经营活动现金流出小计"] + cf_sheet["投资活动现金流出小计"] + cf_sheet["筹资活动现金流出小计"]
        cf_structure_sheet.loc["CFO_OUT"] = cf_sheet["经营活动现金流出小计"] / cf_structure_sheet.loc["现金流出总计"]
        cf_structure_sheet.loc["CFI_OUT"] = cf_sheet["投资活动现金流出小计"] / cf_structure_sheet.loc["现金流出总计"]
        cf_structure_sheet.loc["CFF_OUT"] = cf_sheet["筹资活动现金流出小计"] / cf_structure_sheet.loc["现金流出总计"]
        cf_structure_sheet.to_csv(path_or_buf="./cf_structure_sheet.csv",encoding="gbk")
        return cf_structure_sheet


if __name__ == "__main__":
    antgroup = Quality()
    antgroup.set_period()
    antgroup.result_sheet()
    antgroup.bs_structure_sheet()
    antgroup.is_structure_sheet()
    antgroup.cf_structure_sheet()
